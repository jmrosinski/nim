! Copyright 2016
! United States Government as represented by the National Oceanic and
! Atmospheric Administration.
! No copyright is claimed in the United States under Title 17, U.S.Code.
! All Other Rights Reserved.
!
! Copyright 2016 Regents of the University of Colorado. All rights
! reserved.
! Copyright 2016 Regents of the Colorado State University. All rights
! reserved.

subroutine init 
!*********************************************************************
!
!	Initializes variables and constants
!	Alexander E. MacDonald  11/27/04
!	J. Lee                  September, 2005
!       Jacques Middlecoff      July, 2009 Broke start into start and init.
!       Jacques Middlecoff July 2009 Put variables in calling sequence.
!
!*********************************************************************

use module_control
use module_constants
use module_variables
use globals             ,only: myrank
use module_restart
use permute_io          ,only: read_permute
use ReadNamelist, only: nz

IMPLICIT NONE
#include <gptl.inc>

real              :: zp,zbtm
integer           :: k,ipn,isn,nt,it,ret
integer           :: ios ! return from I/O call
integer,parameter :: unitno=31
character(len=64) :: fn  ! file name

             ret = gptlstart('Init')

! initialize diagnostic array
diagt1d(:,:) = 0.

!....................................................................
!  Runge-Kutta scheme parameters
!
fdt(1)=.5*dt
fdt(2)=.5*dt
fdt(3)=dt

frk(1)=dt/6.
frk(2)=dt/3.
frk(3)=dt/3.
frk(4)=dt/6.

!  numerical diffusion coef

cfh=3.0e3/dt

!$OMP PARALLEL DO PRIVATE(k)
do ipn=ims,ime
  do k=1,nz
    fr (k,ipn) = 0. 
    fur(k,ipn) = 0. 
    fvr(k,ipn) = 0. 
    ftr(k,ipn) = 0. 
  enddo
  do k=0,nz
    fwr(k,ipn) = 0. 
  enddo
enddo

!TBH:  Set initial values of physics quantities that are currently handled 
!TBH:  differently between GRIMS, GFS, and no-physics configurations.  
!TODO:  Clean up this mess.  These values should all be set up during physics 
!TODO:  initialization, possibly further modified if a restart file is read.  
!$OMP PARALLEL DO
do ipn=ips,ipe
  ! TBH:  Set sw2d and lw2d to zero for now to match r1395 baseline 
  ! TBH:  for GRIMs tests.  
  !TODO:  Set sw2d and lw2d to inf.  If GRIMs does not initialize sw2d 
  !TODO:  and lw2d then do not write them from output.F90 when GRIMs 
  !TODO:  is used.  
  sw2d(ipn)      = 0.0
  lw2d(ipn)      = 0.0
  rn2d(ipn)      = 0.0
  rc2d(ipn)      = 0.0
  sn2d(ipn)      = 0.0
enddo

!.......................................................
!Load Define terrain-following vertical coordinate eta
!...................................................................

do k=0,nz
! etae(k)=k*dz
  etae(k)=zh1d(k)
end do

do k=1,nz
  etin(k)=zh1d(k-1)
  etac(k)=.5*(etae(k)+etae(k-1))
end do

spncef=0.
zbtm=2.0e4
do k=0,nz
  zp=.5*pi*(etae(k)-zbtm)/(zt-zbtm)
  if(zp.gt.0.) spncef(k)=(sin(zp))**2
end do

!$OMP PARALLEL DO PRIVATE(k,isn)
do ipn=ips,ipe
  rebb(0,ipn)=rm1d(0)
  tebb(0,ipn)=th1d(0)
  do k=1,nz
    rb(k,ipn)=.5*(rm1d(k-1)+rm1d(k))
    tb(k,ipn)=.5*(th1d(k-1)+th1d(k))
    eb(k,ipn)=((.5*(pm1d(k-1)+pm1d(k)))*1.0e-5)**.286
    rebb(k,ipn)=rm1d(k)
    tebb(k,ipn)=th1d(k)
    trb(k,ipn) = rb(k,ipn)*tb(k,ipn)
  end do
  do isn=1,nprox(ipn)
    do k=1,nz
      reb(k,isn,ipn)=.5*(rm1d(k-1)+rm1d(k))
      teb(k,isn,ipn)=.5*(th1d(k-1)+th1d(k))
    end do
  end do
!	f = coriolis(nip)
!	lat=latitude(nip) degrees
!	lon=longitude(nip) degrees
!	map = map factor at center of edge chords

!  Calculate coriolis acceleration for given icos grid
  f(ipn)=2.*omega*sin(lat(ipn))
  deg_lat(ipn)=raddeg*lat(ipn)
  deg_lon(ipn)=raddeg*lon(ipn)
enddo
!
!....................................................................
!  Load initial state

! read in prognostic variables from "nim_out.dat" to restart model run
! at a particular time step called "ntr"
!
! check if restart is true or not
!
if (RestartBegin .eq. 0) then
  flgrestart = .false.
else
  flgrestart = .true.
endif

!JR Eliminated st2d read since not used in dynamics-only mode

if (.not.flgrestart) then
  fn = 'ini.dat'
  print *,'  Reading ', trim(fn), ' ...'
  open (unitno, file=trim(fn), form="unformatted", action='read', iostat=ios)
  if (ios /= 0) then
    write(6,*) 'init: failure to open input file ', trim(fn), '. stopping'
    stop
  end if
  read (unitno, iostat=ios) it
  if(it.ne.0) stop 'it.ne.0 in init'

  do k=1,nz
    call read_permute (unitno, r(k,:), trim(fn), 'r')
  enddo
  do k=1,nz
    call read_permute (unitno, trp(k,:), trim(fn), 'trp')
  enddo
  do k=1,nz
    call read_permute (unitno, ur(k,:), trim(fn), 'ur')
  enddo
  do k=1,nz
    call read_permute (unitno, vr(k,:), trim(fn), 'vr')
  enddo
  do k=0,nz
    call read_permute (unitno, wr(k,:), trim(fn), 'wr')
  enddo
  do k=1,nz
    call read_permute (unitno, trc1(k,1,:), trim(fn), 'trc1')
  enddo
  close(unitno)
!trp=trp*(1.+0.608*trc1(:,1,:))  ! trp is virtual trp
  do nt=2,ntr1
    trc1(:,nt,:)= 0.
  end do
else ! if (.not.flgrestart) then
!SMS$serial begin
  call get_restart(1,fn)
!SMS$serial end
  print *,'  Reading ', trim(fn), ' ...'
  open (unitno, file=trim(fn), form="unformatted", action='read', iostat=ios)
  if (ios /= 0) then
    write(6,*) 'init: failure to open restart file ', trim(fn), ': stopping'
    stop
  end if
  read (unitno, iostat=ios) RestartTimeUnit
  if (ios /= 0) then
    write(6,*) 'init: failure to read RestartTimeUnit from ', trim(fn)
    call flush(6)
    stop
  end if
  read (unitno, iostat=ios) itatm
  if (ios /= 0) then
    write(6,*) 'init: failure to read itatm from ', trim(fn)
    call flush(6)
    stop
  end if
  write(6,*) '  - Restarting timestep is =',itatm
  do k=1,nz
    call read_permute (unitno, r(k,:), trim(fn), 'r')
  enddo
  do k=1,nz
    call read_permute (unitno, trp(k,:), trim(fn), 'trp')
  enddo
  do k=1,nz
    call read_permute (unitno, ur(k,:), trim(fn), 'ur')
  enddo
  do k=1,nz
    call read_permute (unitno, vr(k,:), trim(fn), 'vr')
  enddo
  do k=0,nz
    call read_permute (unitno, wr(k,:), trim(fn), 'wr')
  enddo
    call read_permute (unitno, rn2d, trim(fn), 'rn2d')
    call read_permute (unitno, rc2d, trim(fn), 'rc2d')
  do nt=1,ntr1
    do k=1,nz
      call read_permute (unitno, trc1(k,:,nt), trim(fn), 'trc1')
    enddo
  enddo
  close(unitno)
  it=itatm
  itsbeg=itatm+1

end if ! (.not.flgrestart) then

! For initial run, perturb the temperature field if this was requested in the namelist  
! NOTE that perturbation is *not* independent of processor count.
if (it == 0 .and. pertlim > 0.d0) then
  call perturb (trp, pertlim)
! Exchange needed to get the perturbed values into halos
  ret = gptlstart('InitExchange')
!sms$exchange (trp)
  ret = gptlstop ('InitExchange')
end if

!$acc update device(ur,vr,wr,trp,r,trb,rb,tb,trc1,tkv,pm1d,spncef)

!This exchange is needed because SMS SERIAL does not do an exchange for a GPU run
             ret = gptlstart('InitExchange')
!sms$exchange (ur,vr,wr,trp,r,trb,rb,tb,trc1)
             ret = gptlstop ('InitExchange')

             ret = gptlstart('Diag')
call diag    (ntr1,ims,ime,ips,ipe,powDouble,p1000,rd,gamma, &
                  kapa,g,qvmin,qwmin,qrmin,area,vol,pm1d,spncef,ca4k, &
                  ca4p,r,rb,trp,trb,tb,trc1,ur,vr,wr, &
                  rp,tp,tr,t,p,qv,qw,qr,u,v, &
                  e,w,tkv)
             ret = gptlstop ('Diag')

!SMS$COMPARE_VAR(tp,"init.F90 - tp")
!SMS$COMPARE_VAR(tr,"init.F90 - tr")
!SMS$COMPARE_VAR(t,"init.F90 - t")
!SMS$COMPARE_VAR(p,"init.F90 - p")
!SMS$COMPARE_VAR(qv,"init.F90 - qv")
!SMS$COMPARE_VAR(qw,"init.F90 - qw")
!SMS$COMPARE_VAR(u,"init.F90 - u")
!SMS$COMPARE_VAR(v,"init.F90 - v")
!SMS$COMPARE_VAR(e,"init.F90 - e")
!SMS$COMPARE_VAR(w,"init.F90 - w")
!SMS$COMPARE_VAR(trc1,"init.F90 - trc1")

if (allocated(globalperm)) then
  ! We are now done with non-decomposed 2D array globalperm
  ! used by read_permute().  
  !TODO:  move globalperm to module permute_io and add init/destroy methods
  deallocate(globalperm)
else
!TODO:  Why does this get printed by every compute task?  
!SMS$IGNORE BEGIN
!  print *,'ERROR task ',myrank,': globalperm not allocated!'
!  call flush(6)
!SMS$IGNORE END
endif

!SMS$INSERT call sms__get_noOverlap(noOverlap)
print*,'noOverlap=',noOverlap

             ret = gptlstop ('Init')
return
end subroutine init

!*********************************************************************
!
!       Nonhydrostatic Icosahedral Model (NIM)
!
!       Design:  Jin-luen Lee and Alexander E. MacDonald (2007-2010)
!       Development Leader:  Jin-luen Lee (2008-2011)
!       Computational Development:  Jacques Middlecoff (2009-2011)
!       GPU Design:  Mark Govett (2009-2011)
!       GPU Development: Mark Govett and Tom Henderson (2009-2011)
!       Documentation:  Alexander E. MacDonald (2010)!
!
!*********************************************************************
