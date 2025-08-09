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

subroutine sponge(ims,ime,ips,ipe,dt,ca4k,ca4p,spncef,r,ur,vr,wr,u,v,w)
!*********************************************************************
!
!	J. Lee                  December, 2005 
!       Jacques Middlecoff May 2010 Put variables in calling sequence.
!	Mark Govett August 2010  GPU parallelization
!
!*********************************************************************
use kinds, only: rt
use ReadNamelist, only: nz
implicit none
integer,intent(IN   ) :: ims,ime,ips,ipe
real(rt),intent(IN   ):: dt
real(rt),intent(IN   ):: ca4k  (1:NZ        )
real(rt),intent(IN   ):: ca4p  (1:NZ        )
real(rt),intent(IN   ):: spncef(0:NZ        ) ! spnge layers coefficient
real(rt),intent(IN   ):: r     (1:NZ,ims:ime)
real(rt),intent(INOUT):: u     (1:NZ,ims:ime)
real(rt),intent(INOUT):: v     (1:NZ,ims:ime)
real(rt),intent(INOUT):: w     (0:NZ,ims:ime)
real(rt),intent(INOUT):: ur    (1:NZ,ims:ime)
real(rt),intent(INOUT):: vr    (1:NZ,ims:ime)
real(rt),intent(INOUT):: wr    (0:NZ,ims:ime)
integer               :: ipn,k,kp1
real(rt)              :: cr

!$acc parallel vector_length(96)
!$OMP PARALLEL DO PRIVATE(k,kp1,cr)
!$acc loop gang
do ipn=ips,ipe
!JR does the kp1 setting disable vectorization? It shouldnt
!$acc loop vector
  do k=1,NZ
    kp1=min(NZ,k+1)
    cr=(ca4k(k)*r(k,ipn)+ca4p(k)*r(kp1,ipn))
    u(k,ipn)= (1._rt-.5_rt*(spncef(k-1)+spncef(k)))*ur(k,ipn)/r(k,ipn)
    v(k,ipn)= (1._rt-.5_rt*(spncef(k-1)+spncef(k)))*vr(k,ipn)/r(k,ipn)
    w(k,ipn)= (1._rt-spncef(k))*wr(k,ipn)/cr
    ur(k,ipn)= r(k,ipn)*u(k,ipn)
    vr(k,ipn)= r(k,ipn)*v(k,ipn)
    wr(k,ipn)= cr*w(k,ipn)
  end do
end do
!$acc end parallel
!$OMP END PARALLEL DO
return
end subroutine sponge

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
