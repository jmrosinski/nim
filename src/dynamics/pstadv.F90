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

subroutine pstadv(ntr1,ims,ime,ips,ipe,npp,nvars,dt,nprox,prox,proxs,vdns,vdnb,sa,area,vol,ca4k,ca4p,sedgvar,trc1)
!*********************************************************************
!
!       Positive definite 3-D finite-volume tracer transport
!	Jin Lee                       Jul, 2011
!
!*********************************************************************
use kinds, only: rt
use ReadNamelist, only: nz
implicit none
integer,intent(IN   ) :: ntr1,ims,ime,ips,ipe,npp,nvars
real(rt),intent(IN   ):: dt
integer,intent(IN   ) :: nprox(          ims:ime  ) ! holds number of proximity points
integer,intent(IN   ) :: prox (      npp,ims:ime  ) ! index of proximity points
integer,intent(IN   ) :: proxs(      npp,ims:ime  ) ! index of proximity sides
real(rt),intent(IN   ):: vdns ( 1:NZ,npp,ims:ime  )
real(rt),intent(IN   ):: vdnb ( 0:NZ,ims:ime,2)
real(rt),intent(IN   ):: sa   (   NZ,npp,ims:ime  )
real(rt),intent(IN   ):: area (      ims:ime      )
real(rt),intent(IN   ):: vol  (   NZ,    ims:ime  ) ! volume for each control volume
real(rt),intent(IN   ):: ca4k (   NZ          )
real(rt),intent(IN   ):: ca4p (   NZ          )
real(rt),intent(IN   ):: sedgvar(1:NZ,npp,nvars,ims:ime)
real(rt),intent(INOUT)::    trc1(1:NZ,ntr1,ims:ime)
real(rt)              :: trfx
real(rt)              :: zfxn ( 0:NZ )
real(rt)              :: zfxp ( 0:NZ )
integer               :: ni  ! tracer index
integer               :: isp,ipp,ipn,k,isn
real(rt)              :: fx1(NZ),fx2(NZ)

!JR For optimal vectorization, should switch the k and isn loops

!$acc parallel num_workers(PAR_WRK) vector_length(VEC_LEN)
!$OMP PARALLEL DO PRIVATE(ni,k,zfxp,zfxn,fx1,fx2,isn,ipp,isp,trfx)
!$acc loop gang private(zfxn,zfxp,fx1,fx2)
do ipn=ips,ipe
!$acc cache(zfxn,zfxp,fx1,fx2)
!$acc loop worker
  do ni=1,ntr1
!$acc loop vector
    do k=1,NZ-1
       zfxp(k)= -.5_rt*(vdnb(k,ipn,2)-abs(vdnb(k,ipn,2)))*(ca4k(k)*trc1(k,ni,ipn)+ca4p(k)*trc1(k+1,ni,ipn))*area(ipn)
       zfxn(k)=  .5_rt*(vdnb(k,ipn,2)+abs(vdnb(k,ipn,2)))*(ca4k(k)*trc1(k,ni,ipn)+ca4p(k)*trc1(k+1,ni,ipn))*area(ipn)
    end do
    zfxn(0) = 0._rt
    zfxp(0) = 0._rt
    zfxn(nz) = 0._rt
    zfxp(nz) = 0._rt
!$acc loop vector
    do k=1,NZ
      fx1(k) = 0._rt
      fx2(k) = 0._rt
    end do
    do isn=1,nprox(ipn)
      ipp = prox(isn,ipn)
      isp = proxs(isn,ipn)
!$acc loop vector
      do k=1,NZ
        trfx = -.5_rt*(vdns(k,isn,ipn)+abs(vdns(k,isn,ipn)))*sedgvar(k,isn,ni,ipn)*sa(k,isn,ipn) &
               +.5_rt*(vdns(k,isp,ipp)+abs(vdns(k,isp,ipp)))*sedgvar(k,isp,ni,ipp)*sa(k,isp,ipp)
        fx1(k) = fx1(k) + max(0._rt, trfx)
        fx2(k) = fx2(k) + max(0._rt,-trfx)
      end do
    end do
!$acc loop vector
    do k=1,NZ
      trc1(k,ni,ipn) = trc1(k,ni,ipn) + ((fx1(k)+zfxp(k)+zfxn(k-1))- &
                                         (fx2(k)+zfxn(k)+zfxp(k-1))) * dt/vol(k,ipn)
!      if(trc1(k,ni,ipn).lt.0) then
!      print *,'negative trc1 (',k,ni,ipn,')=',trc1(k,ni,ipn)
!      stop 'negative trc1 in trc1 '
!      end if
    end do
  end do !ni
enddo !ipn
!$OMP END PARALLEL DO
!$acc end parallel
return
end subroutine pstadv

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
