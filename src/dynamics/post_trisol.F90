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

subroutine post_trisol(nvars,ims,ime,ips,ipe,wr,wrs,ca4k,ca4p,z,nvecb,vdnb,uw8b,bedgvar)
!**************************************************************************
!       Loops prior to call to trisol()
!       Jin Lee                  December, 2008
!**************************************************************************
use kinds, only: rt
use ReadNamelist, only: nz
implicit none
integer,intent(IN   ) :: nvars,ims,ime,ips,ipe
real(rt),intent(IN   ):: wr     (0:NZ,    ims:ime      )
real(rt),intent(IN   ):: wrs    (0:NZ,    ims:ime      )
real(rt),intent(IN   ):: ca4k   (  NZ                  )
real(rt),intent(IN   ):: ca4p   (  NZ                  )
real(rt),intent(IN   ):: z      (0:NZ,    ims:ime      )
real(rt),intent(IN   ):: nvecb  (0:NZ,    ims:ime,    3)
real(rt),intent(  OUT):: vdnb   (0:NZ,    ims:ime,    2)
real(rt),intent(INOUT):: uw8b   (0:NZ,    ims:ime,    3)
real(rt),intent(INOUT):: bedgvar(0:NZ,    ims:ime,nvars)
integer ipn,k,kp1

!$OMP PARALLEL DO PRIVATE(k,kp1)
!$acc parallel vector_length(96)
!$acc loop gang
do ipn=ips,ipe
!JR Does this vectorize? It should.
!$acc loop vector
  do k=0,NZ
    kp1=min(NZ,k+1)
    bedgvar(k,ipn,3)=.5_rt*(wr(k,ipn)+wr(kp1,ipn))
  end do
#ifdef _OPENACC
enddo
!$acc end parallel
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(k,kp1)
!$acc parallel vector_length(96)
!$acc loop gang
do ipn=ips,ipe
#endif
!$acc loop vector
  do k=1,NZ-1
    uw8b(k,ipn,3)=ca4k(k)*bedgvar(k-1,ipn,3)+ca4p(k)*bedgvar(k,ipn,3)
  end do
  uw8b(0 ,ipn,3)= 0._rt

  uw8b(nz,ipn,3)=2._rt*uw8b(nz-1,ipn,3)-uw8b(nz-2,ipn,3)

!$acc loop vector
  do k=0,NZ
    vdnb(k,ipn,1) = 0._rt
    vdnb(k,ipn,2) = uw8b(k,ipn,3)*nvecb(k,ipn,3)
  end do
end do
!$acc end parallel
!$OMP END PARALLEL DO

return
end subroutine post_trisol
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
