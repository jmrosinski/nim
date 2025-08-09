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

subroutine preForce(ims,ime,ips,ipe,npp,nvars,gamma,rd,cfh,nprox,prox,cs,sn,nvecs,sedgvar,sa,e,saoprx,ur,vr,wr,fu0,fv0,fw0)
!*********************************************************************
!Calculates NIM dynamical forcing terms
!
!       Nonhydrostatic Icosahedral Model (NIM)
!  Original subroutine force in QNH:  A.E. MacDonald, Jin Lee -1991
!       Jin Lee                  December, 2008
!  Jacques Middlecoff July 2009 Put K on the inside for the GPU
!  Jacques Middlecoff July 2009 Put variables in calling sequence.
!       Mark Govett August 2010 GPU parallelization
!  Jacques Middlecoff Oct. 2010 Removed SMS
!  Jacques Middlecoff July 2016 Split force into preForce and force.
!*********************************************************************

use kinds, only: rt
use ReadNamelist, only: nz
implicit none
integer ,intent(IN ) :: ims,ime,ips,ipe,npp,nvars
real(rt),intent(IN ) :: gamma,rd,cfh
integer ,intent(IN ) :: nprox  (               ims:ime  ) ! holds number of proximity points
integer ,intent(IN ) :: prox   (     npp,      ims:ime  )
real(rt),intent(IN ) :: cs     (     npp,      ims:ime  )
real(rt),intent(IN ) :: sn     (     npp,      ims:ime  )
real(rt),intent(IN ) :: nvecs  (  NZ,npp,      ims:ime,2)
real(rt),intent(IN ) :: sedgvar(1:NZ,npp,nvars,ims:ime  )
real(rt),intent(IN ) :: sa     (  NZ,npp,      ims:ime  )
real(rt),intent(IN ) :: e      (  NZ,          ims:ime  ) ! Exner func = (p/p0)**kappa
real(rt),intent(IN ) :: saoprx (  NZ,npp,      ims:ime  ) ! diffusion area factor
real(rt),intent(IN ) :: ur     (  NZ,          ims:ime  ) ! u*r, flux variable u times density
real(rt),intent(IN ) :: vr     (  NZ,          ims:ime  ) ! v*r, flux variable v times density
real(rt),intent(IN ) :: wr     (0:NZ,          ims:ime  ) ! w*r, flux variable v times density
real(rt),intent(OUT) :: fu0    (  NZ,          ims:ime  )
real(rt),intent(OUT) :: fv0    (  NZ,          ims:ime  )
real(rt),intent(OUT) :: fw0    (  NZ,          ims:ime  )

integer              :: ipn,isn,k,n,ipp1,ipp2,ipp3,ipp4,ipp5,ippn
real(rt)             :: rhsu(NZ,6),rhsv(NZ,6)
real(rt)             :: sumu(NZ),sumv(NZ),sumw(NZ)

!$OMP PARALLEL DO PRIVATE(k,isn,n,ipp1,ipp2,ipp3,ipp4,ipp5,ippn,rhsu,rhsv, &
!$OMP                     sumu,sumv,sumw) SCHEDULE(runtime)
!$acc parallel num_workers(PAR_WRK) vector_length(VEC_LEN)
!$acc loop gang  private(rhsu,rhsv,sumu,sumv,sumw)
do ipn=ips,ipe
!$acc cache(rhsu,rhsv)
!$acc loop vector
  do k=1,NZ
    fu0(k,ipn) = 0._rt
    fv0(k,ipn) = 0._rt
    fw0(k,ipn) = 0._rt
  enddo
!$acc loop seq
  do isn=1,nprox(ipn)     ! loop thru edges getting fluxes
 !$acc loop vector
    do k=1,NZ
      fu0(k,ipn) = fu0(k,ipn) + sedgvar(k,isn,6,ipn)*nvecs(k,isn,ipn,1)*sa(k,isn,ipn)
      fv0(k,ipn) = fv0(k,ipn) + sedgvar(k,isn,6,ipn)*nvecs(k,isn,ipn,2)*sa(k,isn,ipn)
    end do ! vertical loop
  end do ! end of loop through edges getting fluxes
 !$acc loop vector
  do  k=1,NZ
    fu0(k,ipn) = - gamma*rd*e(k,ipn)*fu0(k,ipn)
    fv0(k,ipn) = - gamma*rd*e(k,ipn)*fv0(k,ipn)
  end do ! vertical loop
  n    = nprox(ipn)
  ipp1 = prox(1,ipn)
  ipp2 = prox(2,ipn)
  ipp3 = prox(3,ipn)
  ipp4 = prox(4,ipn)
  ipp5 = prox(5,ipn)
  ippn = prox(n,ipn)

!DIR$ VECTOR ALWAYS
 !$acc loop vector
  do k=1,NZ
    rhsu(k,1) = cs(1,ipn)*ur(k,ipp1) + sn(1,ipn)*vr(k,ipp1)
    rhsu(k,2) = cs(2,ipn)*ur(k,ipp2) + sn(2,ipn)*vr(k,ipp2)
    rhsu(k,3) = cs(3,ipn)*ur(k,ipp3) + sn(3,ipn)*vr(k,ipp3)
    rhsu(k,4) = cs(4,ipn)*ur(k,ipp4) + sn(4,ipn)*vr(k,ipp4)
    rhsu(k,5) = cs(5,ipn)*ur(k,ipp5) + sn(5,ipn)*vr(k,ipp5)
    rhsu(k,n) = cs(n,ipn)*ur(k,ippn) + sn(n,ipn)*vr(k,ippn)

    rhsv(k,1) =-sn(1,ipn)*ur(k,ipp1) + cs(1,ipn)*vr(k,ipp1)
    rhsv(k,2) =-sn(2,ipn)*ur(k,ipp2) + cs(2,ipn)*vr(k,ipp2)
    rhsv(k,3) =-sn(3,ipn)*ur(k,ipp3) + cs(3,ipn)*vr(k,ipp3)
    rhsv(k,4) =-sn(4,ipn)*ur(k,ipp4) + cs(4,ipn)*vr(k,ipp4)
    rhsv(k,5) =-sn(5,ipn)*ur(k,ipp5) + cs(5,ipn)*vr(k,ipp5)
    rhsv(k,n) =-sn(n,ipn)*ur(k,ippn) + cs(n,ipn)*vr(k,ippn)
  end do

 !$acc loop vector
  do k=1,NZ
    sumu(k) = 0._rt
    sumv(k) = 0._rt
    sumw(k) = 0._rt
  end do

  do isn=1,nprox(ipn)
 !$acc loop vector
    do k=1,NZ
      sumu(k) = sumu(k) + cfh*(rhsu(k,isn)-ur(k,ipn))*saoprx(k,isn,ipn)
      sumv(k) = sumv(k) + cfh*(rhsv(k,isn)-vr(k,ipn))*saoprx(k,isn,ipn)
      sumw(k) = sumw(k) + cfh*(wr(k,prox(isn,ipn))-wr(k,ipn))*saoprx(k,isn,ipn)
    end do
  end do

 !$acc loop vector
  do k=1,NZ
    fu0(k,ipn) = fu0(k,ipn) + sumu(k)
    fv0(k,ipn) = fv0(k,ipn) + sumv(k)
    fw0(k,ipn) = fw0(k,ipn) + sumw(k)
  end do
end do

!$acc end parallel
!$OMP END PARALLEL DO

return
end subroutine preForce

!*********************************************************************
!
!       Design:  Jin-luen Lee and Alexander E. MacDonald (2007-2010)
!       Development Leader:  Jin-luen Lee (2008-2011)
!       Computational Development:  Jacques Middlecoff (2009-2011)
!       GPU Design:  Mark Govett (2009-2011)
!       GPU Development: Mark Govett and Tom Henderson (2009-2011)
!       Documentation:  Alexander E. MacDonald (2010)!
!
!*********************************************************************
