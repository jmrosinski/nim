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

!***************************************************************************
!       Nonhydrostatic Icosahedral Model (NIM)
!
!       Design:  Jin-luen Lee and Alexander E. MacDonald (2007-2010)
!       Development Leader:  Jin-luen Lee (2008-2011)
!       Computational Development:  Jacques Middlecoff (2009-2011)
!       GPU Design:  Mark Govett (2009-2011)
!       GPU Development: Mark Govett and Tom Henderson (2009-2011)
!       Documentation:  Alexander E. MacDonald (2010)!
!       Jacques Middlecoff October 2010 Removed SMS
!
!***************************************************************************
subroutine TimeDiff(ims,ime,ips,ipe,FullStep,fdt, fr,fur,fvr,fwr,ftr, tfr,tfur,tfvr,tfwr,tftr, rs,urs,vrs,wrs,trs, r,ur,vr,wr,trp)
!Time differencing.
use kinds, only: rt
use ReadNamelist, only: nz
implicit none
integer,intent(IN ) :: ims,ime,ips,ipe
logical,intent(IN ) :: FullStep
real(rt),intent(IN ):: fdt                ! partial time step     
real(rt),intent(IN ):: fr  (  NZ,ims:ime)
real(rt),intent(IN ):: fur (  NZ,ims:ime)
real(rt),intent(IN ):: fvr (  NZ,ims:ime)
real(rt),intent(IN ):: fwr (0:NZ,ims:ime)
real(rt),intent(IN ):: ftr (  NZ,ims:ime)
real(rt),intent(IN ):: tfr (  NZ,ims:ime)
real(rt),intent(IN ):: tfur(  NZ,ims:ime)
real(rt),intent(IN ):: tfvr(  NZ,ims:ime)
real(rt),intent(IN ):: tfwr(0:NZ,ims:ime)
real(rt),intent(IN ):: tftr(  NZ,ims:ime)
real(rt),intent(IN )::  rs (  NZ,ims:ime)
real(rt),intent(IN ):: urs (  NZ,ims:ime)
real(rt),intent(IN ):: vrs (  NZ,ims:ime)
real(rt),intent(IN ):: wrs (0:NZ,ims:ime)
real(rt),intent(IN ):: trs (  NZ,ims:ime)
real(rt),intent(OUT)::  r  (  NZ,ims:ime) ! rho (density) prime (kg/m**3) (= r - rb)
real(rt),intent(OUT):: ur  (  NZ,ims:ime) ! u*r, flux variable u times density
real(rt),intent(OUT):: vr  (  NZ,ims:ime) ! v*r, flux variable v times density
real(rt),intent(OUT):: wr  (0:NZ,ims:ime) ! w*r, flux variable w times density
real(rt),intent(OUT):: trp (  NZ,ims:ime) ! flux pot temp prime ( = tr - trb)
integer             :: k,ipn

if(FullStep) then
!$OMP PARALLEL DO PRIVATE(k)
!$acc parallel vector_length(96)
!$acc loop gang
  do ipn=ips,ipe
!$acc loop vector
    do k=1,NZ
       r (k,ipn) =  rs(k,ipn)+fr (k,ipn)
      ur (k,ipn) = urs(k,ipn)+fur(k,ipn)
      vr (k,ipn) = vrs(k,ipn)+fvr(k,ipn)
      wr (k,ipn) = wrs(k,ipn)+fwr(k,ipn)
      trp(k,ipn) = trs(k,ipn)+ftr(k,ipn)
    enddo
    wr(0,ipn) = 0._rt
    wr(NZ,ipn) = 0._rt
  enddo
!$acc end parallel
!$OMP END PARALLEL DO

else

!$OMP PARALLEL DO PRIVATE(k)
!$acc parallel vector_length(96)
!$acc loop gang
  do ipn=ips,ipe
!$acc loop vector
    do k=1,NZ
       r (k,ipn) =  rs(k,ipn)+fdt*tfr (k,ipn)
      ur (k,ipn) = urs(k,ipn)+fdt*tfur(k,ipn)
      vr (k,ipn) = vrs(k,ipn)+fdt*tfvr(k,ipn)
      wr (k,ipn) = wrs(k,ipn)+fdt*tfwr(k,ipn)
      trp(k,ipn) = trs(k,ipn)+fdt*tftr(k,ipn)
    enddo
    wr(0,ipn) = 0._rt
    wr(NZ,ipn) = 0._rt
  enddo
!$acc end parallel
!$OMP END PARALLEL DO
endif
return
end subroutine TimeDiff
