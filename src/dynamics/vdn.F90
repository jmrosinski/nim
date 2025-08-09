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

subroutine vdn(npp,ims,ime,ips,ipe,nprox,prox,proxs,nvecs,nvecb,uw8b, &
               uw8s,vdnb,vdns)
!**********************************************************************
!
!       Calculate the normal wind component for flux calculations
!       Jin Lee                  December, 2008
!       Jacques Middlecoff July 2009 Put variables in calling sequence.
!       Jacques Middlecoff October 2010 Removed SMS
!**********************************************************************
use kinds, only: rt
use ReadNamelist, only: nz
implicit none
integer,intent(IN)  ::         npp,ims,ime,ips,ipe
integer,intent(IN)  :: nprox(         ims:ime)
integer,intent(IN)  :: prox (     npp,ims:ime)
integer,intent(IN)  :: proxs(     npp,ims:ime)
real(rt),intent(IN)  :: nvecs(  NZ,npp,ims:ime,2)
real(rt),intent(IN)  :: nvecb(0:NZ,    ims:ime,3)
real(rt),intent(IN)  :: uw8b (0:NZ,    ims:ime,3)
real(rt),intent(IN)  :: uw8s (1:NZ,npp,2,ims:ime)
real(rt),intent(OUT) :: vdnb (0:NZ,    ims:ime,2)
real(rt),intent(OUT) :: vdns (1:NZ,npp,ims:ime)
integer :: ipn,k,isn

!JR The vdn(0,...) settings are sub-optimal for CPU vectorization

!$acc parallel vector_length(96)
!$acc loop gang
!$OMP PARALLEL DO PRIVATE(isn,k) SCHEDULE(runtime)
do ipn=ips,ipe
  do isn=1,nprox(ipn)
!$acc loop vector
    do k=1,NZ
      vdns(k,isn,ipn) = .5_rt* (uw8s(k,isn,1,ipn)*nvecs(k,isn,ipn,1)  &
                             +  uw8s(k,isn,2,ipn)*nvecs(k,isn,ipn,2)  &
          - uw8s(k,proxs(isn,ipn),1,prox(isn,ipn))*nvecs(k,proxs(isn,ipn),prox(isn,ipn),1)  &
          - uw8s(k,proxs(isn,ipn),2,prox(isn,ipn))*nvecs(k,proxs(isn,ipn),prox(isn,ipn),2)  )
    end do
  end do
  vdnb(0,ipn,1) = 0._rt
  vdnb(0,ipn,2) = uw8b(0,ipn,3)*nvecb(0,ipn,3)
!$acc loop vector
  do k=1,NZ
    vdnb(k,ipn,1) = 0._rt
    vdnb(k,ipn,2) = uw8b(k,ipn,3)*nvecb(k,ipn,3)
  end do
end do  !ipn-loop
!$OMP END PARALLEL DO
!$acc end parallel
return
end subroutine vdn

!**********************************************************************
!
!       Nonhydrostatic Icosahedral Model (NIM)
!
!       Design:  Jin-luen Lee and Alexander E. MacDonald (2007-2010)
!       Development Leader:  Jin-luen Lee (2008-2011)
!       Computational Development:  Jacques Middlecoff (2009-2011)
!       GPU Design:  Mark Govett (2009-2011)
!       GPU Development: Mark Govett and Tom Henderson (2009-2011)
!       Documentation:  Alexander E. MacDonald (2010)
!
!**********************************************************************
