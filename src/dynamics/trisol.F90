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

subroutine trisol(ims,ime,ips,ipe,nvars,cn0,cn1,dt,rd,cv,g,ct4w,caw4tk,caw4tp,ca4p,ca4k,spncef,e,nvecb,tb,rb,st,rp,trp,wr,r,uw8b,vdnb,bedgvar)
!*********************************************************************
!  Jacques Middlecoff SEP 2011 Removed SMS
!  Jacques Middlecoff April 2011 Re-arranged code for the GPU (12%)
!                                Moved ipn loop into trdi.F90 (7%)
!                                Total speedup for CPU serial k32: 20%
!  Jacques Middlecoff April 2011 Ported to the GPU
!  Jacques Middlecoff Nov   2011 Again ported to the GPU
!*********************************************************************

use kinds, only: rt
use ReadNamelist, only: nz
implicit none
integer,intent(IN)    :: ims,ime,ips,ipe,nvars
real(rt),intent(IN)   :: cn0,cn1,dt,rd,cv,g
real(rt),intent(IN)   :: ct4w   (  NZ              )
real(rt),intent(IN)   :: caw4tk (  NZ              )
real(rt),intent(IN)   :: caw4tp (  NZ              )
real(rt),intent(IN)   :: ca4p   (  NZ              )
real(rt),intent(IN)   :: ca4k   (  NZ              )
real(rt),intent(IN)   :: spncef (0:NZ              )
real(rt),intent(IN)   :: e      (  NZ,ims:ime      )
real(rt),intent(IN)   :: nvecb  (0:NZ,ims:ime,3    )
real(rt),intent(IN)   :: rb     (  NZ,ims:ime      )
real(rt),intent(IN)   :: tb     (  NZ,ims:ime      )
real(rt),intent(IN)   :: st     (  NZ,ims:ime      )
real(rt),intent(IN)   :: rp     (  NZ,ims:ime      )
real(rt),intent(INOUT):: trp    (  NZ,ims:ime      )
real(rt),intent(INOUT):: wr     (0:NZ,ims:ime      )
real(rt),intent(INOUT):: r      (  NZ,ims:ime      )
real(rt),intent(OUT)  :: uw8b   (0:NZ,ims:ime,3    )
real(rt),intent(OUT)  :: vdnb   (0:NZ,ims:ime,2    )
real(rt),intent(IN)   :: bedgvar(0:NZ,ims:ime,nvars)

! Local variables

real(rt)              :: w1d    (NZ+1)
real(rt)              :: gama   (NZ  ), alpha
real(rt)              :: aaa    (NZ+1),bbb(NZ+1),ccc(NZ+1),rrr(NZ+1)
real(rt)              :: rtk    (NZ  ),rdk(NZ  )
real(rt)              :: thkp,thkp1,rwk
integer               :: ipn,k,km1,kp1

!$acc parallel num_gangs(10242) vector_length(96) !firstprivate(ct4w,caw4tk,caw4tp,ca4p,ca4k,spncef)
!$acc loop gang private(w1d,aaa,bbb,ccc,rrr,rtk,rdk,gama)  !firstprivate(ct4w,caw4tk,caw4tp,ca4p,ca4k,spncef)
!$OMP PARALLEL DO PRIVATE(ipn,k,kp1,km1,thkp1,thkp,rwk,alpha,w1d,aaa,bbb,ccc,rrr,rtk,rdk,gama)
do ipn=ips,ipe
!$acc cache(w1d,aaa,bbb,ccc,rrr,rdk,rtk,gama)
  uw8b(0 ,ipn,3) = 0._rt
  vdnb(0 ,ipn,1) = 0._rt
!JR Note the following just sets vdnb(0,ipn,2) = 0.
  vdnb(0 ,ipn,2) = uw8b(0,ipn,3)*nvecb(0,ipn,3)
!$acc loop vector
  do k=1,NZ-1
    uw8b(k,ipn,3) = ca4k(k)*bedgvar(k-1,ipn,3)+ca4p(k)*bedgvar(k,ipn,3)
  end do
  uw8b(NZ,ipn,3) = 2._rt*uw8b(NZ-1,ipn,3) - uw8b(NZ-2,ipn,3)

!$acc loop vector
  do k=1,NZ
    vdnb(k,ipn,1) = 0._rt
    vdnb(k,ipn,2) = uw8b(k,ipn,3)*nvecb(k,ipn,3)
  end do
!$acc loop vector
  do k=1,NZ
    rdk(k) = rp(k,ipn) - cn0*ct4w(k)*(vdnb(k,ipn,2) - vdnb(k-1,ipn,2))
    rtk(k) =trp(k,ipn) - cn0*ct4w(k)*(vdnb(k,ipn,2)*bedgvar(k,ipn,5) - vdnb(k-1,ipn,2)*bedgvar(k-1,ipn,5)) &
           + 1._rt*(1._rt - .5_rt*(spncef(k-1) + spncef(k)))*st(k,ipn)*dt
  end do

!JR The NOFUSION directive below is critical--without it ifort 14.0.2 targetting MIC fuses the loop
!JR with one above which results in no vectorization and a 3X slowdown in performance.
!JR This problem only occurs on MIC

!$acc loop vector
!DIR$ NOFUSION
  do k=1,NZ-1  ! k= 0 & nz for BCs
    kp1 = k+1
    km1 = k-1
    thkp1 = .5_rt*( bedgvar(kp1,ipn,6)+bedgvar(k,ipn,6))
    thkp  = .5_rt*( bedgvar(km1,ipn,6)+bedgvar(k,ipn,6))
    rwk   = wr(k,ipn) - cn0*(caw4tk(k)*e(k,ipn) + caw4tp(k)*e(kp1,ipn))*(thkp1 - thkp) 
    rwk   = rwk + cn0*g*((rd/cv)*(thkp1/tb(kp1,ipn) + thkp/tb(k,ipn)) &
            -(ca4k(k)*rp(k,ipn) + ca4p(k)*rp(kp1,ipn)))
    km1 = max(k-1,1)
    aaa(k+1) = -cn1**2*(caw4tk(k)*e(k,ipn) + caw4tp(k)*e(kp1,ipn))*ct4w(k)*bedgvar(k-1,ipn,5)
    bbb(k+1) = 1._rt + cn1**2*(caw4tk(k)*e(k,ipn) + caw4tp(k)*e(kp1,ipn))*(ct4w(k) + ct4w(kp1))*bedgvar(k,ipn,5)
    ccc(k+1) = -cn1**2*(caw4tk(k)*e(k,ipn) + caw4tp(k)*e(kp1,ipn))*ct4w(kp1)*bedgvar(k+1,ipn,5)
    rrr(k+1) = rwk - cn1*(caw4tk(k)*e(k,ipn) + caw4tp(k)*e(kp1,ipn))*(rtk(kp1)-rtk(k)) 
    aaa(k+1) = aaa(k+1) - cn1**2*ct4w(k)*g*ca4k(k)*((rd/cv)/tb(k,ipn)*bedgvar(k-1,ipn,5)-1._rt)
    bbb(k+1) = bbb(k+1) + cn1**2*(-(g*rd/cv)*(ca4p(k)*ct4w(kp1)/tb(kp1,ipn) - &
                                              ca4k(k)*ct4w(k  )/tb(k,  ipn))*bedgvar(k,ipn,5))
    bbb(k+1) = bbb(k+1) - cn1**2*g*(ca4k(k)*ct4w(k)- ca4p(k)*ct4w(kp1))
    ccc(k+1) = ccc(k+1) - cn1**2*ct4w(kp1)*g*ca4p(k)*(-(rd/cv)/tb(kp1,ipn)*bedgvar(k+1,ipn,5) + 1._rt)
    rrr(k+1) = rrr(k+1) + cn1*g*((rd/cv)*(ca4p(k)/tb(kp1,ipn)*rtk(kp1) + ca4k(k)/tb(k,ipn)*rtk(k)) &
               -(ca4k(k)*rdk(k) + ca4p(k)*rdk(kp1)))
  end do

  aaa(1) = 0._rt
  bbb(1) = 1._rt
  ccc(1) = 0._rt
  rrr(1) = 0._rt
  aaa(NZ+1) = 0._rt
  bbb(NZ+1) = 1._rt
  ccc(NZ+1) = 0._rt
  rrr(NZ+1) = 0._rt
!  call trdslv0(nz,ims,ime,ips,ipe,aaa,bbb,ccc,rrr,w1d,gama)
!     OBTAIN THE LU-DECOMPOSITION
  w1d(1)  = rrr(1)/bbb(1)
  gama(1) = ccc(1)/bbb(1)
!$acc loop seq
  DO k=2,NZ
    alpha = 1._rt/(bbb(k)-aaa(k)*gama(k-1))
    gama(k) = ccc(k)*alpha
    w1d(k) = (rrr(k)-aaa(k)*w1d(k-1))*alpha
  ENDDO

!     SOLVE

  w1d(NZ+1) = 0._rt
!$acc loop seq
  DO k=NZ,1,-1
    w1d(k) = w1d(k) - gama(k)*w1d(k+1)
  ENDDO
  wr(0,ipn) = w1d(1)
!$acc loop vector
  do k=1,NZ
    wr(k,ipn) = w1d(k+1)
  end do
!$acc loop vector
  do k=1,NZ
    trp(k,ipn) = rtk(k) - cn1*ct4w(k)*(wr(k,ipn)*bedgvar(k,ipn,5) - wr(k-1,ipn)*bedgvar(k-1,ipn,5))
    r(k,ipn)   = rb(k,ipn) + rdk(k) - cn1*ct4w(k)*(wr(k,ipn)-wr(k-1,ipn))
  end do
end do  !ipn-loop
!$OMP END PARALLEL DO
!$acc end parallel

return
end subroutine trisol

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
