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

#define G3DCONTAINSPRX

subroutine input
!*********************************************************************
!
!	Reads the initial variables and constants.
!	Alexander E. MacDonald  11/27/04
!	J. Lee                  September, 2005
!       Jacques Middlecoff      July, 2009 Separated start into start, input and init.
!
!*********************************************************************
  use module_control
  use module_constants
  use module_variables
  use permute_io      ,only: read_permute
  use amtx            ,only: compute_amtx
  use kinds           ,only: rt, sp
  use globals         ,only: myrank,parallelBuild
  use ReadNamelist, only: nz

  IMPLICIT NONE
#include <gptl.inc>

!TODO: Allocate unit number with units.F90
  integer, parameter :: unitno = 28
  integer :: nprocs=1
  integer :: ipn,k,ret,isn,i,ia,ierr
  integer :: tmpi4(nip)   ! integer array read from disk
  real(sp):: ca4kr4(nz),ca4pr4(nz),caw4tkr4(nz),caw4tpr4(nz),ct4wr4(nz)

#ifdef _OPENMP
  integer, external :: omp_get_max_threads
#endif
  CHARACTER(80) :: FileName

  ret = gptlstart ('Input')

!..................................................................
!Load Icosahedral grid values

!!SMS$SERIAL (<zs,zsc,zsb,zbx,zbxb,OUT> : default=ignore)  BEGIN
!FileName="topo.dat"
!open(unit=unitno, file=TRIM(FileName), form="unformatted", status='old', action='read', err=70)
!read(unitno, err=90) zs,zsc,zsb,zbx,zbxb
!close(unitno)
!!SMS$SERIAL END

! Per Jin set these all to zero for aqua-planet
!JR Changed to zero only the interior region.
!SMS$PARALLEL (dh,ipn) BEGIN
!$OMP PARALLEL DO
  do ipn=1,nip
    zs(ipn)     = 0.
    zsc(:,ipn)  = 0.
    zsb(:,ipn)  = 0.
    zbx(ipn)    = 0.
    zbxb(:,ipn) = 0.
  end do

  FileName = "g3d.dat"
  open(unit=unitno, file=TRIM(FileName), form="unformatted", status='old', action='read', iostat=ierr)
  if (ierr.ne.0) then
    write (*,'(a,a)') 'input: error opening ',trim(filename)
    stop
  endif
  print *,'  Reading ',TRIM(FileName),' ...'
!TODO:  optimize SERIAL INOUT
  call read_permute (unitno, lat, FileName, 'lat')
  call read_permute (unitno, lon, FileName, 'lon')
  do k=0,nz
    call read_permute (unitno, z(k,:), FileName, 'z')
  enddo

  do k=1,nz
    call read_permute (unitno, zm(k,:), FileName, 'zm')
  enddo

!SMS$SERIAL (default=ignore) BEGIN
  do isn=1,npp
    read(unitno, iostat=ierr) tmpi4     ! prox already set in start.F90
    if (ierr.ne.0) then
      write (*,'(a,i0)') 'input: error reading prox.',isn
      stop
    endif
  enddo
  read(unitno, iostat=ierr) tmpi4     ! nprox already set in start.F90
  if (ierr.ne.0) then
    write (*,'(a)') 'input: error reading nprox'
    stop
  endif

  do isn=1,npp
    read(unitno, iostat=ierr) tmpi4     ! proxs already set in start.F90
    if (ierr.ne.0) then
      write (*,'(a,i0)') 'input: error reading proxs',isn
      stop
    endif
  enddo
!SMS$SERIAL END
#ifdef G3DCONTAINSPRX
!JR NOTE: order of isn,i loops DIFFERS from datagen/g3d.F90
!JR But, serial NIM works this way so the order must be correct.
!JR Reversing the order causes NIM to fail.
  do isn=1,npp
    do i=1,2
      call read_permute (unitno, prox_xy(isn,i,:), FileName, 'prox_xy')
    enddo
  enddo
#endif

  call read_permute (unitno, area, FileName, 'area')

  do isn=1,npp
    call read_permute (unitno, cs(isn,:), FileName, 'cs')
  enddo

  do isn=1,npp
    call read_permute (unitno, sn(isn,:), FileName, 'sn')
  enddo

  do isn=1,npp
    do i=1,2
      call read_permute (unitno, xyc(i,isn,:), FileName, 'xyc')
    enddo
  enddo

  do isn=1,npp
    do i=1,2
      call read_permute (unitno, xys(i,isn,:), FileName, 'xys')
    enddo
  enddo

  do isn=1,npp
    do i=1,2
      call read_permute (unitno, xyb(i,isn,:), FileName, 'xyb')
    enddo
  enddo

  do isn=1,npp
    do k=0,nz
      call read_permute (unitno, zc(k,isn,:), FileName, 'zc')
    enddo
  enddo

  do isn=1,npp
    do k=1,nz
      call read_permute (unitno, zcs(k,isn,:), FileName, 'zcs')
    enddo
  enddo

  do i=1,2
    do isn=1,npp
      do k=1,nz
        call read_permute (unitno, nvecs(k,isn,:,i), FileName, 'nvecs')
      enddo
    enddo
  enddo

  do i=1,3
    do k=0,nz
      call read_permute (unitno, nvecb(k,:,i), FileName, 'nvecb')
    enddo
  enddo

  do isn=1,npp
    do k=1,nz
      call read_permute (unitno, sa(k,isn,:), FileName, 'sa')
    enddo
  enddo

  do k=1,nz
    call read_permute (unitno, vol(k,:), FileName, 'vol')
  enddo

  do isn=1,npp
    do k=1,nz
      call read_permute (unitno, saoprx(k,isn,:), FileName, 'saoprx')
    enddo
  enddo

!SMS$SERIAL (<ca4k,ca4p,caw4tk,caw4tp,ct4w,OUT> : default=ignore)  BEGIN
  write(6,*)'reading ca4k,etc.'
  read(unitno, iostat=ierr) ca4kr4,ca4pr4,caw4tkr4,caw4tpr4,ct4wr4
  if (ierr.ne.0) then
    write (*,'(a)') 'input: error reading ca4k,etc.'
    stop
  endif
  ca4k   = ca4kr4
  ca4p   = ca4pr4
  caw4tk = caw4tkr4
  caw4tp = caw4tpr4
  ct4w   = ct4wr4
!SMS$SERIAL END
  close(unitno)

!$acc enter data create(vdns,vdnb,urs,vrs,wrs,trs,rs,bedgvar,efr,tefr,tfr,tfur,tfvr,tfwr,tftr,fr,fur,fvr,fwr,ftr,fu0,fv0,fw0,qv,qw,qr,rp,tp,tr,t,p,u,v,e,w, &
!$acc                   nvecs,nvecb,sa,proxs,cs,sn,vol,nprox,lat,lon,area,ca4k,ca4p,z,zm,ur,vr,wr,trp,r,trb,rb,tb,trc1,tkv,pm1d,spncef,                     &
!$acc                   ct4w,caw4tk,caw4tp,f,amtx1,amtx2,uw8s,uw8b,zc,zcs,zb3d,reb,rebb,teb,tebb,xys,xyb,xyc,prox,saoprx,div,sedgvar,st,su,sv,sw)

!$acc update device(nvecs,nvecb,sa,proxs,cs,sn,vol,nprox,lat,lon,area,ca4k,ca4p,z,zm)

! Update halos of constant arrays.  Cannot rely on automatic halo update from 
! file read due to prox not yet being set up.
! For the GPU everything exchanged here must be copied out after the exchange if the halo is used on the CPU.
! 11/21/10: Only sa and nprox need to be exchanged (halos are used).
!           The other 8 variables are exchanged for possible future use of halos.
!           For GPU runs: None of the ten variables have their halos used on the CPU.
!                         All 10 variables are copied to the CPU for possible future use of their halos on the CPU.
  ret = gptlstart('InitExchange')
!sms$exchange (nvecs,nvecb,sa,nprox,cs,sn,lat,lon,area,vol,z,zm)
  ret = gptlstop ('InitExchange') 
!$acc update host(nvecs,nvecb,z,zm,sa,proxs,nprox,cs,sn,lat,lon)

  if (read_amtx) then
    FileName="amx.dat"
    open(unitno, file=TRIM(FileName), form="unformatted", status='old', action='read', iostat=ierr)
    if (ierr.ne.0) then
      write (*,'(a,a)') 'input: error opening ',trim(filename)
      stop
    endif

    print *,'  Reading amtx arrays from file=',trim(FileName),' ...'

    ret = gptlstart ('read_amtx')
    do k=1,nz
      do ia=1,(nob+1)*nbf
        call read_permute (unitno, amtx1(k,ia,:), FileName, 'amtx1')
      enddo
    enddo

    do k=1,nz-1 
      do ia=1,(nob+1)*nbf
        call read_permute (unitno, amtx2(k,ia,:), FileName, 'amtx2')
      enddo
    enddo
    ret = gptlstop ('read_amtx')

    close(unitno)

  else
    write(6,*) 'Computing amtx arrays...'
    ret = gptlstart ('compute_amtx')
    call compute_amtx ()
    ret = gptlstop ('compute_amtx')
  end if
!SMS$PARALLEL END

!TBH:  First, test with first-touch loops only, then add this call
!TBH:  and replace divides with multiplies in sovlveiThLS*() routines.  
#ifdef SOLVEITHLS_RECIPROCAL
  call invert_A_diag(nz,nob,ims,ime,ips,ipe,amtx1,amtx2)
#endif

!...................................................................
!Load mean profile

  FileName = "STD.txt"
!SMS$SERIAL (<zh1d,tm1d,pm1d,rm1d,th1d,OUT> : default=ignore)  BEGIN
  open (unit=5, file=FileName, form="formatted", status='old', action='read', iostat=ierr)
  if (ierr.ne.0) then
    write (*,'(a,a)') 'input: error opening ',trim(filename)
    stop
  endif
  print *,'  Reading ',TRIM(FileName),' ...'
  do k=0,nz
    read (5,*,iostat=ierr) zh1d(k),tm1d(k),pm1d(k),rm1d(k)
    if (ierr.ne.0) then
      write (*,'(a)') 'input: error reading zh1d(k),tm1d(k),pm1d(k),rm1d(k)'
      stop
    endif
    th1d(k) = tm1d(k)*(p1000/pm1d(k))**kapa
  end do ! k-loop
  close(5)
!SMS$SERIAL END

!...................................................................
!  Print input variables

!sms$comm_size(nprocs)
  print"('NIM Global Model')"
  print *, " "
  call datetime
  print *, " "
  print"(' Grid Level                           ',I9              )",glvl
  print"(' Number of Processors:                ',I9,' processors')",nprocs
#ifdef _OPENMP
  print"(' omp_get_max_threads():               ',I9,' threads')",omp_get_max_threads()
#endif
  print"(' Global size:                         ',I9,' points'    )",nip
!SMS$insert print"(' Halo size:                           ',I9,' points'    )",haloSize
!           print"(' Forecast duration:                   ',I9,' days, ',I0,' hours')",numday,numhour
  print"(' Vertical resolution:                 ',I9,' levels'    )",nz
  print"(' Length of time step:                 ',F9.2,' seconds' )",dt
  print"(' Number of time steps:                ',I9,' timesteps' )",nts
  print"(' Output interval:                     ',I9,' timesteps')",ArchvStep
  print"(' Physics package:                     ',A9              )",trim(physics)
  print"(' writeOutput:                         ',L9              )",writeOutput
  print"(' powDouble:                           ',L9              )",powDouble
  print"(' vdmint_combine:                      ',I9              )",vdmint_combine
  print"(' preExchangeBarrier                   ',L9              )",preExchangeBarrier
  print"(' postExchangeBarrier                  ',L9              )",postExchangeBarrier
  print"(' TimeInitExchanges                    ',L9              )",TimeInitExchanges
  print"(' pin_to_single_core                   ',L9              )",pin_to_single_core
  print"(' root_on_socket1                      ',L9              )",root_on_socket1
  print"(' outputBarrier                        ',L9              )",outputBarrier
  print"(' zeroSMStimers                        ',L9              )",zeroSMStimers
  print"(' read_amtx                            ',L9              )",read_amtx
  print *,' '
  print *,' '

  ret = gptlstop ('Input')
  return
end subroutine input

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
