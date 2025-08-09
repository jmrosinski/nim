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

subroutine vdmfinish(npp,nvars,ims,ime,ips,ipe,nprox,rebb,tebb,reb,ca4k,ca4p,bedgvar,sedgvar,uw8s,uw8b)
!**************************************************************************
!       Final loops of vdm.F90 after calls to vdmints()
!       Jin Lee                  December, 2008
!**************************************************************************
use kinds, only: rt
use ReadNamelist, only: nz
implicit none
integer,intent(IN   ) :: npp,nvars,ims,ime,ips,ipe
integer,intent(IN   ) :: nprox  (         ims:ime      )
real(rt),intent(IN   ) :: rebb   (0:nz,    ims:ime      )
real(rt),intent(IN   ) :: tebb   (0:nz,    ims:ime      )
real(rt),intent(IN   ) :: reb    (  nz,npp,ims:ime      )
real(rt),intent(IN   ) :: ca4k   (  nz                  )
real(rt),intent(IN   ) :: ca4p   (  nz                  )
real(rt),intent(INOUT) :: bedgvar(0:nz,    ims:ime,nvars)
real(rt),intent(INOUT) :: sedgvar(  nz,npp,nvars,ims:ime)
real(rt),intent(  OUT) :: uw8s   (  nz,npp,2,ims:ime)
real(rt),intent(  OUT) :: uw8b   (0:nz,    ims:ime,    3)
integer ipn,k,isn

!DIR$ ASSUME_ALIGNED reb:64, ca4k:64, ca4p:64, sedgvar:64
!$acc parallel vector_length(96)!
!$acc loop gang
!$OMP PARALLEL DO PRIVATE(k,isn) SCHEDULE(runtime)
do ipn=ips,ipe
  bedgvar(0,ipn,3) = 0._rt
  bedgvar(0,ipn,4) = rebb(0,ipn) + bedgvar(0,ipn,4)
  bedgvar(0,ipn,5) = tebb(0,ipn) + bedgvar(0,ipn,5)
!$acc loop vector
  do k=1,nz
    bedgvar(k,ipn,4) = max(-.95_rt*rebb(k,ipn), bedgvar(k,ipn,4))
    bedgvar(k,ipn,4) = rebb(k,ipn) + bedgvar(k,ipn,4)
    bedgvar(k,ipn,5) = tebb(k,ipn) + bedgvar(k,ipn,5)
  end do
  do isn=1,nprox(ipn)
!$acc loop vector
    do k=1,nz
      sedgvar(k,isn,4,ipn) = reb    (k,isn,ipn  )+sedgvar(k,isn,4,ipn)
      uw8s   (k,isn,1,ipn) = sedgvar(k,isn,4,ipn)*sedgvar(k,isn,1,ipn)
      uw8s   (k,isn,2,ipn) = sedgvar(k,isn,4,ipn)*sedgvar(k,isn,2,ipn) 
   end do
  end do
  uw8b(0,ipn,1) = 0._rt
  uw8b(0,ipn,2) = 0._rt
  uw8b(0,ipn,3) = 0._rt
!$acc loop vector
  do k=1,nz
    uw8b(k,ipn,1) = 0._rt
    uw8b(k,ipn,2) = 0._rt
  end do
!$acc loop vector
  do k=1,nz-1
    uw8b(k,ipn,3) = bedgvar(k,ipn,4)*(ca4k(k)*bedgvar(k-1,ipn,3) + ca4p(k)*bedgvar(k,ipn,3))
  end do
  uw8b(nz,ipn,3) = 0._rt
end do
!$OMP END PARALLEL DO
!$acc end parallel

return
end subroutine vdmfinish
!**************************************************************************
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
!**************************************************************************
