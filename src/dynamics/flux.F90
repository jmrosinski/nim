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

subroutine flux(npp,ims,ime,ips,ipe,nvars,nprox,prox,proxs,vdns,vdnb,sa,area,vol,ca4k,ca4p,sedgvar,bedgvar,tfur,tfvr,tfwr,tftr,tefr)
!**********************************************************************
!Calculate finite-volume 3-D integrated quantities.
!	Jin Lee            Dec 2008
!       Jacques Middlecoff May 2010  Put variables in calling sequence.
!	Mark Govett        Aug 2010  GPU parallelization
!       Jacques Middlecoff Oct 2010  Combined the five calls to flux into one flux routine.
!       Jacques Middlecoff Oct 2010  Removed SMS
!**********************************************************************

use kinds, only: rt
use ReadNamelist, only: nz
implicit none
integer ,intent(IN   ):: npp,ims,ime,ips,ipe,nvars
integer ,intent(IN   ):: nprox  (         ims:ime  )!holds number of proximity points
integer ,intent(IN   ):: prox   (     npp,ims:ime  )!index of proximity points
integer ,intent(IN   ):: proxs  (     npp,ims:ime  )!index of proximity sides
real(rt),intent(IN   ):: vdns   (1:NZ,npp,ims:ime  )
real(rt),intent(IN   ):: vdnb   (0:NZ,    ims:ime,2)
real(rt),intent(IN   ):: sa     (  NZ,npp,ims:ime  )
real(rt),intent(IN   ):: area   (         ims:ime  )
real(rt),intent(IN   ):: vol    (  NZ,    ims:ime  )!volume for each control volume
real(rt),intent(IN   ):: ca4k   (  NZ              )
real(rt),intent(IN   ):: ca4p   (  NZ              )
real(rt),intent(IN   ):: sedgvar(1:NZ,npp,nvars,ims:ime)
real(rt),intent(IN   ):: bedgvar(0:NZ,    ims:ime,nvars)
real(rt),intent(INOUT):: tfur   (1:NZ    ,ims:ime      )
real(rt),intent(INOUT):: tfvr   (1:NZ    ,ims:ime      )
real(rt),intent(INOUT):: tfwr   (0:NZ    ,ims:ime      )
real(rt),intent(INOUT):: tftr   (1:NZ    ,ims:ime      )
real(rt),intent(  OUT):: tefr   (  NZ,npp,ims:ime      )
real(rt) :: vnk,vnkm1
integer  :: isp,ipp,ipn,isn,k,kp1
real(rt) :: fx(NZ)
real(rt) :: fx1(NZ)
real(rt) :: fx2(NZ)
real(rt) :: fx5(NZ)
real(rt) :: fipn,fipp,upfx1,upfx2,upfx3

!$OMP PARALLEL DO PRIVATE(k,isn,ipp,isp,fx1,fx2,fx5,fx,kp1,vnkm1,vnk,upfx1,upfx2,upfx3) SCHEDULE(runtime)
!$acc parallel num_workers(PAR_WRK) vector_length(VEC_LEN)
!$acc loop gang private(fx,fx1,fx2,fx5)
do ipn=ips,ipe
!$acc cache(fx,fx1,fx2,fx5)
!$acc loop vector
  do k=1,NZ
    fx1(k) = 0._rt
    fx2(k) = 0._rt
    fx5(k) = 0._rt
  end do
!$acc loop seq
  do isn=1,nprox(ipn)     ! loop thru edges getting fluxes
    ipp=prox( isn,ipn)
    isp=proxs(isn,ipn)
!$acc loop vector
    do k=1,NZ
!Compared to the old calculation of fx1,fx2 and fx5 below,
!this new calculation has digits accuracy of 5,5,3,5,6 for -O3 or -O0
!with physics commented out.
!For the old and new code alike, -O3 compared to -O0 is bitwise exact!
!     fipn = .5_rt*(vdns(k,isn,ipn)+abs(vdns(k,isn,ipn)))*sa(k,isn,ipn)
!     fipp = .5_rt*(vdns(k,isp,ipp)+abs(vdns(k,isp,ipp)))*sa(k,isp,ipp)
!     fx1  = fx1 + fipn*sedgvar(k,isn,1,ipn) - fipp*sedgvar(k,isp,1,ipp)
!     fx2  = fx2 + fipn*sedgvar(k,isn,2,ipn) - fipp*sedgvar(k,isp,2,ipp)
!     fx5  = fx5 + fipn*sedgvar(k,isn,5,ipn) - fipp*sedgvar(k,isp,5,ipp)
      tefr (k,isn,ipn) = .5_rt*(vdns(k,isp,ipp)+abs(vdns(k,isp,ipp)))  &
                       - .5_rt*(vdns(k,isn,ipn)+abs(vdns(k,isn,ipn)))
      fx1(k) = fx1(k)  + .5_rt*(vdns(k,isn,ipn)+abs(vdns(k,isn,ipn)))*sedgvar(k,isn,1,ipn)*sa(k,isn,ipn) &
                       - .5_rt*(vdns(k,isp,ipp)+abs(vdns(k,isp,ipp)))*sedgvar(k,isp,1,ipp)*sa(k,isp,ipp)
      fx2(k) = fx2(k)  + .5_rt*(vdns(k,isn,ipn)+abs(vdns(k,isn,ipn)))*sedgvar(k,isn,2,ipn)*sa(k,isn,ipn) &
                       - .5_rt*(vdns(k,isp,ipp)+abs(vdns(k,isp,ipp)))*sedgvar(k,isp,2,ipp)*sa(k,isp,ipp)
      fx5(k) = fx5(k)  + .5_rt*(vdns(k,isn,ipn)+abs(vdns(k,isn,ipn)))*sedgvar(k,isn,5,ipn)*sa(k,isn,ipn) &
                       - .5_rt*(vdns(k,isp,ipp)+abs(vdns(k,isp,ipp)))*sedgvar(k,isp,5,ipp)*sa(k,isp,ipp)
    end do
  end do ! end of loop through edges getting fluxes

!$acc loop vector
  do k=1,NZ
    tfur(k,ipn) = tfur(k,ipn) - fx1(k)/vol(k,ipn)
    tfvr(k,ipn) = tfvr(k,ipn) - fx2(k)/vol(k,ipn)
    tftr(k,ipn) = tftr(k,ipn) - fx5(k)/vol(k,ipn)
  end do

!$acc loop vector
  do k=1,NZ-1
    fx(k) = 0._rt
  end do

!!$acc loop seq
  do isn=1,nprox(ipn)     ! loop thru edges getting fluxes
    ipp=prox(isn,ipn)
    isp=proxs(isn,ipn)
!$acc loop vector
    do k=1,NZ-1
      fx(k) = fx(k) +.5_rt*(ca4k(k)*(vdns(k,isn,ipn)+abs(vdns(k,isn,ipn)))+ca4p(k)*(vdns(k+1,isn,ipn)+abs(vdns(k+1,isn,ipn)))) &
             *sedgvar(k,isn,3,ipn)*(ca4p(k)*sa(k,isn,ipn)+ca4p(k)*sa(k+1,isn,ipn))   &
                    -.5_rt*(ca4k(k)*(vdns(k,isp,ipp)+abs(vdns(k,isp,ipp)))+ca4p(k)*(vdns(k+1,isp,ipp)+abs(vdns(k+1,isp,ipp)))) & 
             *sedgvar(k,isp,3,ipp)*(ca4k(k)*sa(k,isp,ipp)+ca4p(k)*sa(k+1,isp,ipp))
    end do
  end do ! end of loop through edges getting fluxes
!$acc loop vector
  do k=1,NZ-1
    tfwr(k,ipn) = tfwr(k,ipn) - fx(k)/vol(k,ipn)
  end do
  tfwr(NZ,ipn) = 2._rt*tfwr(NZ-1,ipn) - tfwr(NZ-2,ipn)
!JR Split k loop into 2 pieces because previous kp1=min(k+1,NZ) wouldnt vectorize
!$acc loop vector
  do k=1,NZ-1
    kp1   = k+1
    vnkm1 =.5_rt*(vdnb(k  ,ipn,2)+vdnb(k-1,ipn,2))
    vnk   =.5_rt*(vdnb(kp1,ipn,2)+vdnb(k  ,ipn,2))
    upfx1 = vdnb(k,ipn,2)*bedgvar(k,ipn,1)-vdnb(k-1,ipn,2)*bedgvar(k-1,ipn,1)
    upfx2 = vdnb(k,ipn,2)*bedgvar(k,ipn,2)-vdnb(k-1,ipn,2)*bedgvar(k-1,ipn,2)
    upfx3 = vnk          *bedgvar(k,ipn,3)-vnkm1          *bedgvar(k-1,ipn,3)
    tfur(k,ipn) = tfur(k,ipn) - upfx1*area(ipn)/vol(k,ipn)
    tfvr(k,ipn) = tfvr(k,ipn) - upfx2*area(ipn)/vol(k,ipn)
    tfwr(k,ipn) = tfwr(k,ipn) - upfx3*area(ipn)/(ca4k(k)*vol(k,ipn)+ca4p(k)*vol(kp1,ipn))
  end do
  k = nz
  kp1 = nz
  vnkm1 =.5_rt*(vdnb(k  ,ipn,2)+vdnb(k-1,ipn,2))
  vnk   =.5_rt*(vdnb(kp1,ipn,2)+vdnb(k  ,ipn,2))
  upfx1 = vdnb(k,ipn,2)*bedgvar(k,ipn,1)-vdnb(k-1,ipn,2)*bedgvar(k-1,ipn,1)
  upfx2 = vdnb(k,ipn,2)*bedgvar(k,ipn,2)-vdnb(k-1,ipn,2)*bedgvar(k-1,ipn,2)
  upfx3 = vnk          *bedgvar(k,ipn,3)-vnkm1          *bedgvar(k-1,ipn,3)
  tfur(k,ipn) = tfur(k,ipn) - upfx1*area(ipn)/vol(k,ipn)
  tfvr(k,ipn) = tfvr(k,ipn) - upfx2*area(ipn)/vol(k,ipn)
  tfwr(k,ipn) = tfwr(k,ipn) - upfx3*area(ipn)/(ca4k(k)*vol(k,ipn)+ca4p(k)*vol(kp1,ipn))
end do ! ipn-loop
!$acc end parallel
!$OMP END PARALLEL DO

return
end subroutine flux

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
