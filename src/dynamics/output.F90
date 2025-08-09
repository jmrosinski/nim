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

subroutine output(nip,itsout,minmaxPrintInterval,itsend,&
                  numvars,var_list,r,wr,trp,u,v,w,t,tp,p,rp,ur,vr)
!*********************************************************************
!
!       Output program for non-hydrostatic global model
!       Modified from sgm_output.F90 by Alexander E. MacDonald  12/27/2004
!       J. Lee                  December, 2008
!       Jacques Middlecoff July 2009 Put variables in calling sequence.
!       Jacques Middlecoff Sept 2009 Converted to 1 file per variable
!
!*********************************************************************

  use icosio,             only: icosio_end_frame,icosio_out
  use module_restart          ! rad_inthr
  use module_constants,   only: pi
  use module_control,     only: dt,filenamelen,ips,ipe,ims,ime,OutputByHour,physics, &
                                ntr1,ArchvStep,step_per_unit,ArchvTimeUnit
  use module_header,      only: header
  use module_variables,   only: e, su, sv, qr, qw, qv, tr, trc1, global_inv_perm
  use module_variables,   only: st2d,rn2d,rc2d,sw2d,lw2d,sfc_shflx,sfc_lhflx
  use module_its2time
  use globalutils,        only: mininfo, maxinfo
  use ReadNamelist,       only: nz,outputBarrier
  use kinds, only: rt, sp
  implicit none

  integer,intent (IN) :: itsout,minmaxPrintInterval,itsend,nip,numvars
  character(128),intent (IN) :: var_list

!SMS$DISTRIBUTE(dh,1) BEGIN
  real(rt)             :: tmp2d   (nip)
  real(sp)             :: tmp2dinv(nip)
!SMS$DISTRIBUTE END

!SMS$DISTRIBUTE(dh,2) BEGIN
  real(rt),intent (IN) :: r  (  nz,nip) ! density (kg/m**3) r for "rho"
  real(rt),intent (IN) :: wr (0:nz,nip) ! w*r, flux variable w times density
  real(rt),intent (IN) :: trp(  nz,nip) ! flux pot temp prime ( = tr - trb)	
  real(rt),intent (IN) :: u  (  nz,nip) ! zonal wind (m/s)
  real(rt),intent (IN) :: v  (  nz,nip) ! meridional wind (m/s)
  real(rt),intent (IN) :: w  (0:nz,nip) ! vertical velocity (m/s)
  real(rt),intent (IN) :: t  (  nz,nip) ! potential temperature (K)
  real(rt),intent (IN) :: tp (  nz,nip) ! potential temperature prime (k) (= t - tb)
  real(rt),intent (IN) :: p  (0:nz,nip) ! pressure (pascals)
  real(rt),intent (IN) :: rp (  nz,nip) ! rho (density) prime (kg/m**3) (= r - rb)
  real(rt),intent (IN) :: ur (  nz,nip) ! u*r, flux variable u times density
  real(rt),intent (IN) :: vr (  nz,nip) ! v*r, flux variable v times density
  real(rt)             :: var(nz,nip),var1(0:nz,nip)
!SMS$DISTRIBUTE END
#include <gptl.inc>

  integer             :: ipn,k,nvi,time,ret,time_hr
  character(len=filenamelen),external :: filename
  character(len=8)                    :: iname
  character(len=64)                   :: oname
  CHARACTER(len=4)    :: var_names(200), var_tmp_nm
  INTEGER             :: idx, i, minmaxPrintTimeStep

if(outputBarrier) then
!sms$barrier
endif
ret = gptlstart ('Output')
  
  if(mod(itsout,ArchvStep)==0.or.itsout==itsend) then

!$acc update host (trp,ur,vr,wr,r,t,p,qr,qw,qv,u,v,w)

  time_hr = (itsout+.5_rt)*dt/3600
  time    = nint(its2time(itsout))

    if (itsout /= 0 .and. mod(dble(itsout*dt),3600._rt*rad_inthr)==0._rt) then
!SMS$SERIAL BEGIN
      iname = 'rst_atm_'
      call get_fname(iname,8,time,oname)
      open(unit=29,file=oname,form="unformatted", action='write')
      write(29)   ArchvTimeUnit
      write(29)   itsout
!SMS$SERIAL END

      do k=1,nz
!TODO:  need a new write_permute() subroutine
        do ipn=ips,ipe
          tmp2d(ipn)=r(k,ipn)
        enddo
!TODO:  Need a "NONE" clause for SERIAL to use with tmp2dinv since it does *not* 
!TODO:  need to be scattered or gathered!  
!SMS$SERIAL (<tmp2d,tmp2dinv,IN> : default=ignore) BEGIN
        do ipn=1,nip
          tmp2dinv(ipn)=tmp2d(global_inv_perm(ipn))
        enddo
        write(29) tmp2dinv
!SMS$SERIAL END
      enddo

      do k=1,nz
        do ipn=ips,ipe
          tmp2d(ipn)=trp(k,ipn)
        enddo
!SMS$SERIAL (<tmp2d,tmp2dinv,IN> : default=ignore) BEGIN
        do ipn=1,nip
          tmp2dinv(ipn)=tmp2d(global_inv_perm(ipn))
        enddo
        write(29) tmp2dinv
!SMS$SERIAL END
      enddo

      do k=1,nz
        do ipn=ips,ipe
          tmp2d(ipn)=ur(k,ipn)
        enddo
!SMS$SERIAL (<tmp2d,tmp2dinv,IN> : default=ignore) BEGIN
        do ipn=1,nip
          tmp2dinv(ipn)=tmp2d(global_inv_perm(ipn))
        enddo
        write(29) tmp2dinv
!SMS$SERIAL END
      enddo

      do k=1,nz
        do ipn=ips,ipe
          tmp2d(ipn)=vr(k,ipn)
        enddo
!SMS$SERIAL (<tmp2d,tmp2dinv,IN> : default=ignore) BEGIN
        do ipn=1,nip
          tmp2dinv(ipn)=tmp2d(global_inv_perm(ipn))
        enddo
        write(29) tmp2dinv
!SMS$SERIAL END
      enddo

      do k=0,nz
        do ipn=ips,ipe
          tmp2d(ipn)=wr(k,ipn)
        enddo
!SMS$SERIAL (<tmp2d,tmp2dinv,IN> : default=ignore) BEGIN
        do ipn=1,nip
          tmp2dinv(ipn)=tmp2d(global_inv_perm(ipn))
        enddo
        write(29) tmp2dinv
!SMS$SERIAL END
      enddo

      do ipn=ips,ipe
        tmp2d(ipn)=rn2d(ipn)
      enddo
!SMS$SERIAL (<tmp2d,tmp2dinv,IN> : default=ignore) BEGIN
      do ipn=1,nip
        tmp2dinv(ipn)=tmp2d(global_inv_perm(ipn))
      enddo
      write(29) tmp2dinv
!SMS$SERIAL END

      do ipn=ips,ipe
        tmp2d(ipn)=rc2d(ipn)
      enddo
!SMS$SERIAL (<tmp2d,tmp2dinv,IN> : default=ignore) BEGIN
      do ipn=1,nip
        tmp2dinv(ipn)=tmp2d(global_inv_perm(ipn))
      enddo
      write(29) tmp2dinv
!SMS$SERIAL END

      do nvi=1,ntr1
        do k=1,nz
          do ipn=ips,ipe
            tmp2d(ipn)=trc1(k,nvi,ipn)
          enddo
!SMS$SERIAL (<tmp2d,tmp2dinv,IN> : default=ignore) BEGIN
          do ipn=1,nip
            tmp2dinv(ipn)=tmp2d(global_inv_perm(ipn))
          enddo
          write(29) tmp2dinv
!SMS$SERIAL END
        enddo
      enddo
    endif

      idx = 0
      var_tmp_nm(:) = achar(32)

do i=1,numvars
   var_tmp_nm(1:4) = var_list(idx+1:idx+4)
   var_names(i) = var_tmp_nm
   idx = idx + 5
   if (var_names(i).eq.'pZZZ') then
      var1 = p
   elseif (var_names(i).eq.'wZZZ') then
      var1 = w
   elseif (var_names(i).eq.'wrZZ') then
      var1 = wr
   elseif (var_names(i).eq.'urZZ') then
      var = ur
   elseif (var_names(i).eq.'uZZZ') then
      var = u
   elseif (var_names(i).eq.'vrZZ') then
      var = vr
   elseif (var_names(i).eq.'vZZZ') then
      var = v
   elseif (var_names(i).eq.'trpZ') then
      var = trp
   elseif (var_names(i).eq.'trZZ') then
      var = tr 
   elseif (var_names(i).eq.'tpZZ') then
      var = tp
   elseif (var_names(i).eq.'tZZZ') then
      var = t
   elseif (var_names(i).eq.'rpZZ') then
      var = rp
   elseif (var_names(i).eq.'qvZZ') then
      var = qv
   elseif (var_names(i).eq.'qwZZ') then
      var = qw
   elseif (var_names(i).eq.'qrZZ') then
      var = qr
   elseif (var_names(i).eq.'rZZZ') then
      var = r
   else
      print*, var_names(i), ' is undefined inside the var_list in NIMnamelist.'
      print*, 'Please check with output.F90 if the variable is being defined.'
      stop
   end if
   if (var_names(i).eq.'pZZZ'.or.var_names(i).eq.'wZZZ'.or.var_names(i).eq.'wrZZ')then
      call icosio_out(itsout,time,var_names(i),var1,filename(var_names(i),time),header(var_names(i),nz+1,itsout,time_hr,ArchvTimeUnit))
   else
      call icosio_out(itsout,time,var_names(i),var ,filename(var_names(i),time),header(var_names(i),nz  ,itsout,time_hr,ArchvTimeUnit))
   end if
end do
    call icosio_end_frame(itsout)
  end if

  minmaxPrintTimeStep = minmaxPrintInterval * step_per_unit

  if (mod(itsout,minmaxPrintTimeStep) == 0) then
!$acc update host (trp,ur,vr,wr,r,t,p,qr,qw,qv,u,v,w)

    ret = gptlstart ('minmaxprint')
    call mininfo (r,   'r  ')
    call maxinfo (r,   'r  ')

    call mininfo (qv,  'qv ')
    call maxinfo (qv,  'qv ')

    call mininfo (qw,  'qw ')
    call maxinfo (qw,  'qw ')

    call mininfo (qr,  'qr ')
    call maxinfo (qr,  'qr ')

    call mininfo (u,   'u  ')
    call maxinfo (u,   'u  ')

    call mininfo (v,   'v  ')
    call maxinfo (v,   'v  ')

    call mininfo (w,   'w  ')
    call maxinfo (w,   'w  ')

    call mininfo (trp, 'trp')
    call maxinfo (trp, 'trp')
    ret = gptlstop ('minmaxprint')
  end if

  ret = gptlstop ('Output')
  return
end subroutine output

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
