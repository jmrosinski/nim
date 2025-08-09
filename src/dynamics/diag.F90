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

subroutine diag(ntr, ims, ime, ips, ipe, powDouble, p1000, rd, gamma, &
                   kapa, g, qvmin, qwmin, qrmin, area, vol, pm1d, spncef, ca4k, &
                   ca4p, r, rb, trp, trb, tb, trc1, ur, vr, wr, &
                   rp, tp, tr, t, p, qv, qw, qr, u, v, &
                   e, w, tkv)
!*********************************************************************
!
!       Calculate diagnostic and decouple variables
!       Jin-Luen Lee December, 2008
!       Jacques Middlecoff July 2009 Put variables in calling sequence.
!       Mark Govett        July 2010 GPU parallelization
!
!*********************************************************************
  use kinds, only: rt
  use ReadNamelist, only: nz
  implicit none

  integer ,intent (IN)    :: ntr,ims,ime,ips,ipe
  logical ,intent (IN)    :: powDouble          ! do ** in double
  real(rt),intent (IN)    :: p1000,rd,gamma,kapa,g,qvmin,qwmin,qrmin
  real(rt),intent (IN)    :: area(    ims:ime) ! the area of cell polygon
  real(rt),intent (IN)    :: vol(  NZ,ims:ime) ! volume for each control volume
  real(rt),intent (IN)    :: ca4k  (1:NZ)       
  real(rt),intent (IN)    :: ca4p  (1:NZ)       
  real(rt),intent (IN)    :: spncef(0:NZ)        
  real(rt),intent (IN)    :: pm1d(0:NZ)        ! mean pressure
  real(rt),intent (IN)    :: r  (  NZ,ims:ime) ! density (kg/m**3) r for "rho"
  real(rt),intent (IN)    :: rb (  NZ,ims:ime) ! density basic state (kg/m**3)
  real(rt),intent (IN)    :: trp(  NZ,ims:ime) ! flux pot temp prime ( = tr - trb)	
  real(rt),intent (IN)    :: trb(  NZ,ims:ime) ! flux potential temperature basic state
  real(rt),intent (IN)    :: tb (  NZ,ims:ime) ! potential temperature basic state (k)
  real(rt),intent (INOUT) :: trc1( NZ,ntr,ims:ime) ! tracer coupled variables (e.g., qv times density)
  real(rt),intent (IN)    :: ur (  NZ,ims:ime) ! u*r, flux variable u times density
  real(rt),intent (IN)    :: vr (  NZ,ims:ime) ! v*r, flux variable v times density
  real(rt),intent (IN)    :: wr (0:NZ,ims:ime) ! w*r, flux variable w times density
  real(rt),intent (OUT)   :: rp (  NZ,ims:ime) ! rho (density) prime (kg/m**3) (= r - rb)
  real(rt),intent (OUT)   :: tp (  NZ,ims:ime) ! potential temperature prime (k) (= t - tb)
  real(rt),intent (OUT)   :: tr (  NZ,ims:ime) ! t*r, flux variable t times density
  real(rt),intent (OUT)   :: t  (  NZ,ims:ime) ! potential temperature (K)
  real(rt),intent (OUT)   :: p  (0:NZ,ims:ime) ! pressure (pascals)
  real(rt),intent (OUT)   :: qv (  NZ,ims:ime) ! specific humidity for water vapor
  real(rt),intent (OUT)   :: qw (  NZ,ims:ime) ! specific humidity for cloud water
  real(rt),intent (OUT)   :: qr (  NZ,ims:ime) ! specific humidity for cloud water
  real(rt),intent (OUT)   :: u  (  NZ,ims:ime) ! zonal wind (m/s)
  real(rt),intent (OUT)   :: v  (  NZ,ims:ime) ! meridional wind (m/s)
  real(rt),intent (OUT)   :: e  (  NZ,ims:ime) ! Exner func = (p/p0)**kappa
  real(rt),intent (OUT)   :: tkv(  NZ,ims:ime) ! virtual temperature 
  real(rt),intent (OUT)   :: w  (0:NZ,ims:ime) ! vertical velocity (m/s)
  real(rt)                   :: pdel( NZ) ! Change in pressure (pascals)
  integer  :: ipn,k
  real(rt) :: eh
  real(rt) :: qwk
  real(rt) :: term(NZ)

!$acc parallel num_workers(PAR_WRK) vector_length(VEC_LEN)
!$acc loop gang private(pdel,term)
!$OMP PARALLEL DO PRIVATE(k,qwk,eh,pdel,term) SCHEDULE (runtime)
  do ipn=ips,ipe
!$acc cache(pdel,term)
!$acc loop vector
    do k=1,NZ
      eh=.5_rt*(spncef(k-1)+spncef(k))
      rp(k,ipn)=  r (k,ipn)- rb(k,ipn) 
      qwk = .5_rt* ( ( abs(max(-0.1_rt, abs(rp(k,ipn)/rb(k,ipn)) ) ) - max(-0.1_rt, rp(k,ipn)/rb(k,ipn) ) ) & 
                 -( abs(min( 0.1_rt, abs(rp(k,ipn)/rb(k,ipn)) ) ) + min( 0.1_rt, rp(k,ipn)/rb(k,ipn) ) )   )
      tp(k,ipn)=(trp(k,ipn)- rp(k,ipn)*tb(k,ipn))*( (1._rt-eh)/r(k,ipn) + eh*(1._rt+qwk)/rb(k,ipn) )
      t (k,ipn)=  tp(k,ipn)+ tb(k,ipn)  ! t is air potentail virtual temperature (t=tr/r)
      tr(k,ipn)= trp(k,ipn)+trb(k,ipn)  ! tr is tv*rho
      u (k,ipn)=  ur(k,ipn)*( (1._rt-eh)/r(k,ipn) + eh*(1._rt+qwk)/rb(k,ipn) )
      v (k,ipn)=  vr(k,ipn)*( (1._rt-eh)/r(k,ipn) + eh*(1._rt+qwk)/rb(k,ipn) )
      qv(k,ipn)= max( trc1(k,1,ipn)/r(k,ipn), qvmin)
      trc1(k,1,ipn)= r(k,ipn)*qv(k,ipn)
      pdel(k) = g*r(k,ipn)*(1._rt+qv(k,ipn))*vol(k,ipn)/area(ipn)
    end do
    p(NZ,ipn)=pm1d(NZ)
!$acc loop seq
    do k=NZ-1,0,-1
      p(k,ipn)=p(k+1,ipn)+pdel(k+1)
    end do

!JR This can become a Fortran "if" once the custom intrinsics are actually
!JR available.
#ifdef CUSTOM_INTRINSICS
!$acc loop vector
      do k=1,NZ
        e(k,ipn) = slowpowf(p1000*slowpowf (rd*tr(k,ipn)*1.e-5_rt, gamma)*1.e-5_rt, kapa)
        eh       = slowpowf(0.5_rt*(p(k-1,ipn)+p(k,ipn))/p1000,kapa)    ! exner fn for hydrostatic comp
        tkv(k,ipn) = t(k,ipn)*eh
      end do
#else
!JR Rewrote loops to isolate "**" from other computations. Problem is, host cant vectorize
!JR "**" with "-fp-model precise" compiler flag. -fp-model precise is required for bitwise
!JR reproducibility across differing MPI task counts. Why KNC gives reproducible results without
!JR this flag is a mystery. The new construct speeds up host on this routine by around 30%, with 
!JR little to no effect on MIC. Though host is still much slower than MIC, presumably due to no
!JR vectorization of "**". Total speedup of dynamics in symmetric mode is around 3%. 
!JR "NOFUSION" is required because otherwise ifort wants to fuse the loop with a previous one,
!JR which doesnt vectorize due to a dependency!

!$acc loop vector
!DIR$ NOFUSION
      do k=1,NZ
        term(k) = rd*tr(k,ipn)*1.e-5_rt
      end do
!$acc loop vector
      do k=1,NZ
        if (powDouble) then
          term(k) = dble(term(k))**dble(gamma)
        else
          term(k) = term(k)**gamma
        end if
      end do
!$acc loop vector
      do k=1,NZ
        term(k) = (p1000*term(k))*1.e-5_rt
      end do
!$acc loop vector
      do k=1,NZ
        if (powDouble) then
          e(k,ipn) = dble(term(k))**dble(kapa)
        else
          e(k,ipn) = term(k)**kapa
        end if
      end do

!$acc loop vector
      do k=1,NZ
        term(k) = 0.5_rt*(p(k-1,ipn)+p(k,ipn))/p1000
      end do
!$acc loop vector
      do k=1,NZ
        if (powDouble) then
          term(k) = dble(term(k))**dble(kapa)    ! exner fn for hydrostatic comp
        else
          term(k) = term(k)**kapa                ! exner fn for hydrostatic comp
        end if
      end do
!$acc loop vector
      do k=1,NZ
        tkv(k,ipn) = t(k,ipn)*term(k)
      end do
#endif

    if (ntr >= 3) then
!$acc loop vector
      do k=1,NZ
        !This IF test replaces MAX to prevent negative zero on the GPU
        qwk = trc1(k,2,ipn)/r(k,ipn)
        if(qwk <= qwmin) then
          qw(k,ipn)= qwmin
        else
          qw(k,ipn)= qwk
        endif
        trc1(k,2,ipn)= r(k,ipn)*qw(k,ipn)
        qr(k,ipn)    = max( trc1(k,3,ipn)/r(k,ipn), qrmin)
        trc1(k,3,ipn)= r(k,ipn)*qr(k,ipn)
      end do
    else
!$acc loop vector
      do k=1,NZ
        qw(k,ipn)=0._rt
        qr(k,ipn)=0._rt
      end do
    end if

    w(0,ipn)= 0._rt
!$acc loop vector
    do k=1,NZ-1
      w(k,ipn) = wr(k,ipn)/(ca4k(k)*r(k,ipn) + ca4p(k)*r(k+1,ipn))
    end do
    w(NZ,ipn)=0._rt
  end do ! ipn-loop
!$OMP END PARALLEL DO
!$acc end parallel

  return
end subroutine diag
