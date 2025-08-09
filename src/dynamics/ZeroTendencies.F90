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

subroutine Zerotendencies(ims,ime,ips,ipe,fr,fur,fvr,fwr,ftr)
use kinds, only: rt
use ReadNamelist, only: nz
implicit none
integer,intent (IN ) ::          ims,ime,ips,ipe
real(rt),intent (OUT):: fr (  NZ,ims:ime)  ! forcing of density (r)
real(rt),intent (OUT):: fur(  NZ,ims:ime)  ! forcing of ur
real(rt),intent (OUT):: fvr(  NZ,ims:ime)  ! forcing of vr
real(rt),intent (OUT):: fwr(0:NZ,ims:ime)  ! forcing of wr
real(rt),intent (OUT):: ftr(  NZ,ims:ime)  ! forcing of tr
integer              :: k,ipn

!$acc parallel vector_length(96)
!$OMP PARALLEL DO PRIVATE (k)
!$acc loop gang
do ipn=ips,ipe
  fwr(0,ipn) = 0._rt
!$acc loop vector
  do k=1,NZ
    fur(k,ipn) = 0._rt
    fvr(k,ipn) = 0._rt
    ftr(k,ipn) = 0._rt
    fr (k,ipn) = 0._rt
    fwr(k,ipn) = 0._rt
  end do
end do
!$acc end parallel
!$OMP END PARALLEL DO

return

end subroutine Zerotendencies

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
