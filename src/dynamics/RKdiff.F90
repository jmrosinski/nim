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
!       Runge-Kutta differencing.
!       Note that important optimizations from Jacques Middlecoff committed 
!       with r863 are no longer present.  From Jacques commit message:  
!  r863 | jacques | 2011-12-14 20:45:57 -0700 (Wed, 14 Dec 2011) | 6 lines
!  Reorganized and optimized RKdiff reducing the local storage by 12X.
!  On the CPU RKdiff is 1.9X faster.
!  Jacques May 2013:
!     Removed cr since cr=1 always and demoted tefrp,efrp,efrn,tefrn to scalars.
!     Put everything in one REGION statement since Mark fixed the problem with IF statements.
!***************************************************************************
subroutine RKdiff(npp,ims,ime,ips,ipe,FullStep,frk,dt,rhomin,nprox,sa,vol,tefr,efr,tfr,tfur,tfvr,tfwr,tftr,fr,fur,fvr,fwr,ftr)
use kinds, only: rt
use ReadNamelist, only: nz
implicit none
integer,intent(IN   ) :: npp,ims,ime,ips,ipe
logical,intent(IN   ) :: FullStep
real(rt),intent(IN   ):: frk,dt,rhomin
integer,intent(IN   ) :: nprox(         ims:ime) ! holds number of proximity points
real(rt),intent(IN   ):: sa   (  NZ,npp,ims:ime)
real(rt),intent(IN   ):: vol  (  NZ,    ims:ime) ! volume for each control volume
real(rt),intent(IN   ):: tefr (  NZ,npp,ims:ime)
real(rt),intent(IN   ):: tfur (  NZ,    ims:ime) ! tendency function of ur
real(rt),intent(IN   ):: tfvr (  NZ,    ims:ime) ! tendency function of vr
real(rt),intent(IN   ):: tfwr (0:NZ,    ims:ime) ! tendency function of wr
real(rt),intent(IN   ):: tftr (  NZ,    ims:ime) ! tendency function of tr
real(rt),intent(INOUT):: efr  (  NZ,npp,ims:ime) ! edge forcing of density (r)
real(rt),intent(INOUT):: tfr  (  NZ,    ims:ime) ! tendency function of density (r)
real(rt),intent(INOUT):: fr   (  NZ,    ims:ime) ! forcing of density (r)
real(rt),intent(INOUT):: fur  (  NZ,    ims:ime) ! forcing of ur
real(rt),intent(INOUT):: fvr  (  NZ,    ims:ime) ! forcing of vr
real(rt),intent(INOUT):: fwr  (0:NZ,    ims:ime) ! forcing of wr
real(rt),intent(INOUT):: ftr  (  NZ,    ims:ime) ! forcing of tr
real(rt)              :: fx(NZ)
real(rt)              :: tefrp,efrp,efrn,tefrn
integer               :: ipn,k,isn

!$acc parallel vector_length(96)
!$OMP PARALLEL DO PRIVATE(k,isn,fx,tefrn,tefrp,efrn,efrp) SCHEDULE(runtime)
!$acc loop gang private(fx)
do ipn=ips,ipe
!$acc cache(fx)
  if(FullStep) then
    fwr(0,ipn) = frk*tfwr(0,ipn) 
!$acc loop vector
    do k=1,NZ
      fur(k,ipn) = frk*tfur(k,ipn) 
      fvr(k,ipn) = frk*tfvr(k,ipn) 
      fwr(k,ipn) = frk*tfwr(k,ipn) 
      ftr(k,ipn) = frk*tftr(k,ipn)
    end do
    do isn=1,nprox(ipn)
!$acc loop vector
      do k=1,NZ
        efr(k,isn,ipn) = frk*tefr(k,isn,ipn)
      end do
    enddo
  else
    fwr(0,ipn) = fwr(0,ipn) + frk*tfwr(0,ipn) 
!$acc loop vector
    do k=1,NZ
      fur(k,ipn) = fur(k,ipn) + frk*tfur(k,ipn) 
      fvr(k,ipn) = fvr(k,ipn) + frk*tfvr(k,ipn) 
      fwr(k,ipn) = fwr(k,ipn) + frk*tfwr(k,ipn) 
      ftr(k,ipn) = ftr(k,ipn) + frk*tftr(k,ipn)
    end do
    do isn=1,nprox(ipn)
!$acc loop vector
      do k=1,NZ
        efr(k,isn,ipn) = efr(k,isn,ipn) + frk*tefr(k,isn,ipn)
      end do
    enddo
  endif

!$acc loop vector
  do k=1,NZ
    fx(k) = 0._rt
  end do
  do isn=1,nprox(ipn)     ! loop thru edges getting fluxes
!$acc loop vector
    do k=1,NZ
      tefrn = max(0._rt,-tefr (k,isn,ipn))
      tefrp = max(0._rt, tefr (k,isn,ipn))
      fx(k) = fx(k) + (tefrp-tefrn)*sa(k,isn,ipn)
    end do
  end do
!$acc loop vector
  do k=1,NZ
    tfr(k,ipn) = tfr(k,ipn) + fx(k)/vol(k,ipn)
  end do

  if( .not. FullStep .and. (frk < 0.2_rt*dt) ) then
!$acc loop vector
    do k=1,NZ
      fx(k) = 0._rt
    end do
    do isn=1,nprox(ipn)     ! loop thru edges getting fluxes
!$acc loop vector
      do k=1,NZ
        efrn = max(0._rt,-efr (k,isn,ipn))
        efrp = max(0._rt, efr (k,isn,ipn))
        fx(k) = fx(k) + (efrp-efrn)*sa(k,isn,ipn)
      end do
    end do
!$acc loop vector
    do k=1,NZ
      fr(k,ipn) = fr(k,ipn) + fx(k)/vol(k,ipn)
    end do
  end if
end do
!$acc end parallel 
!$OMP END PARALLEL DO

return
end subroutine RKdiff
