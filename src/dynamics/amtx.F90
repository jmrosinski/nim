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

!*********************************************************************
!  J. Lee                  January, 2009
!    Define three-dimensional control volume boxes          
!  J. Rosinski             November, 2013
!    Moved from offline program into NIM model
!*********************************************************************

module amtx
!JR The following commented out use associations will work when newppp is implemented. 
!JR For now the use associations need to be in the contained subroutines
!JR nz is actually in ReadNamelist but module_control uses it
!  use module_control,   only: nip, nob, nbf, nz
!  use module_variables, only: amtx1, amtx2, lat, lon, prox, nprox, proxs, z, zm, prox_xy, &
!                              ca4k, ca4p
  use module_control, only: nob, nbf
  use infnan, only: inf

  implicit none

  private
  public :: compute_amtx

  type Site3d
    real :: x
    real :: y
    real :: z
  end type Site3d

  integer, parameter :: nBasisFuns = nbf
  integer, parameter :: nobs = nob
  integer, parameter :: noa = (nob+1)*nbf

contains

  subroutine compute_amtx ()
    use ReadNamelist, only: nz
    use module_control,   only: nip, nob, ips, ipe
    use module_variables, only: amtx1, amtx2
    TYPE(Site3d) :: stf(nobs)

!JR Moved "as" into compute_amtx to ease OMP-ization of ipn loops
    real :: as(noa)
    integer :: k, ipn

!JR Unfortunately something below preSolveiThLS prevents OMP from working
!JR when applied to the following ipn loops. So turn off. Probably the
!JR code inside of nim_lapack.f that says "if(first)then" because the
!JR occasional failure is at nim_lapack.f:2815 with a divide by zero.
    do k=1,nz-1
      do ipn=ips,ipe
        CALL assignObsLoc(k,k+1,ipn,stf)
        CALL constructiThObsMat(stf, 1, as)
        CALL preSolveiThLS(1, as)
        amtx1(k,1:noa,ipn) = as(1:noa)
      end do
    end do

    do ipn=ips,ipe
      CALL assignObsLoc(nz-1,nz,ipn,stf)
      CALL constructiThObsMat(stf, 1, as)
      CALL preSolveiThLS(1, as)
      amtx1(nz,1:noa,ipn) = as(1:noa)
    end do

    do k=1,nz-1
      do ipn=ips,ipe
        CALL assignObsLow(k,k+1,ipn,stf)
        CALL constructiThObsMat(stf, 1, as)
        CALL preSolveiThLS(1, as)
        amtx2(k,1:noa,ipn) = as(1:noa)
      end do
    end do

    return
  end subroutine compute_amtx

!**********************************************************************
!
!       Assign stencil locations in box (k,ipn)
!       Jin Lee                  February, 2009
!
!**********************************************************************
!
  subroutine assignObsLoc(k1,k2,ipn,stf)
    use module_control,   only: nob
    use module_variables, only: prox, nprox, zm, prox_xy, &
                                ca4k, ca4p
    integer, intent(IN) :: k1
    integer, intent(IN) :: k2
    integer, intent(IN) :: ipn 
    TYPE(Site3d), intent(inout) :: stf(nob)

    real :: zbase

    zbase = zm(k1,ipn)

    stf(1)%x = prox_xy(1,1,ipn) 
    stf(1)%y = prox_xy(1,2,ipn) 
    stf(1)%z = (zm(k1,prox(1,ipn))-zbase)

    stf(2)%x = prox_xy(2,1,ipn) 
    stf(2)%y = prox_xy(2,2,ipn) 
    stf(2)%z = (zm(k1,prox(2,ipn))-zbase)

    stf(3)%x = prox_xy(3,1,ipn) 
    stf(3)%y = prox_xy(3,2,ipn) 
    stf(3)%z = (zm(k1,prox(3,ipn))-zbase)

    stf(4)%x = prox_xy(5,1,ipn) 
    stf(4)%y = prox_xy(5,2,ipn) 
    stf(4)%z = (zm(k1,prox(5,ipn))-zbase)

    stf(5)%x = prox_xy(2,1,ipn) 
    stf(5)%y = prox_xy(2,2,ipn) 
    stf(5)%z = (zm(k1+1,prox(2,ipn))-zbase)

    stf(6)%x = prox_xy(3,1,ipn) 
    stf(6)%y = prox_xy(3,2,ipn) 
    stf(6)%z = (zm(k1+1,prox(3,ipn))-zbase)

    stf(7)%x = prox_xy(4,1,ipn) 
    stf(7)%y = prox_xy(4,2,ipn)
    stf(7)%z = (zm(k1+1,prox(4,ipn))-zbase)
    
    stf(8)%x = prox_xy(nprox(ipn),1,ipn) 
    stf(8)%y = prox_xy(nprox(ipn),2,ipn)
    stf(8)%z = (zm(k1+1,prox(nprox(ipn),ipn))-zbase)

    stf(9)%x = 0.
    stf(9)%y = 0.
    stf(9)%z = ((ca4k(k1)*zm(k1,ipn) + ca4p(k1)*zm(k1+1,ipn)) - zbase)

    return
  end subroutine assignObsLoc

!**********************************************************************
!
!       Assign stencil locations in box (k,ipn)
!       Jin Lee                  February, 2009
!
!**********************************************************************
!
  subroutine assignObsLow(k1,k2,ipn,stf)
    use module_control,   only: nob, nbf
    use module_variables, only: prox, nprox, z, prox_xy
    integer, intent(IN) :: k1
    integer, intent(IN) :: k2
    integer, intent(IN) :: ipn 
    TYPE(Site3d), intent(inout) :: stf(nob)

    real :: zbase
!
!    zbase=.5*(z(k1,ipn)+z(k2,ipn))
    zbase = z(k1,ipn)

    stf(1)%x = prox_xy(1,1,ipn)
    stf(1)%y = prox_xy(1,2,ipn) 
    stf(1)%z = (z(k1,prox(1,ipn))-zbase)

    stf(2)%x = prox_xy(2,1,ipn) 
    stf(2)%y = prox_xy(2,2,ipn) 
    stf(2)%z = (z(k1,prox(2,ipn))-zbase)

    stf(3)%x = prox_xy(3,1,ipn) 
    stf(3)%y = prox_xy(3,2,ipn) 
    stf(3)%z = (z(k1,prox(3,ipn))-zbase)

    stf(4)%x = prox_xy(5,1,ipn) 
    stf(4)%y = prox_xy(5,2,ipn) 
    stf(4)%z = (z(k1,prox(5,ipn))-zbase)

    stf(5)%x = prox_xy(2,1,ipn) 
    stf(5)%y = prox_xy(2,2,ipn) 
    stf(5)%z = (z(k1+1,prox(2,ipn))-zbase)

    stf(6)%x = prox_xy(3,1,ipn) 
    stf(6)%y = prox_xy(3,2,ipn) 
    stf(6)%z = (z(k1+1,prox(3,ipn))-zbase)

    stf(7)%x = prox_xy(4,1,ipn) 
    stf(7)%y = prox_xy(4,2,ipn) 
    stf(7)%z = (z(k1+1,prox(4,ipn))-zbase)

    stf(8)%x = prox_xy(nprox(ipn),1,ipn) 
    stf(8)%y = prox_xy(nprox(ipn),2,ipn) 
    stf(8)%z = (z(k1+1,prox(nprox(ipn),ipn))-zbase)

    stf(9)%x = 0.
    stf(9)%y = 0.
    stf(9)%z = (.5*(z(k1,ipn)+z(k2,ipn))-zbase)

    return
  end subroutine assignObsLow

!---------------------------------------------------------------------------
!  This file contains the classes that implement the least square fitting. 
!  using quadratic polynomial functions
!
!  
!  Ning Wang, Nov 13 2013
!  
!---------------------------------------------------------------------------

  SUBROUTINE constructiThObsMat(st3d, ith, as)
    TYPE(Site3d), intent(in) :: st3d(nObs)
    INTEGER, intent(in) :: ith
    real, intent(inout) :: as(noa)

    INTEGER :: i, j, ithPtStart, sizeA

    sizeA = nObs * nBasisFuns
    ithPtStart = (ith - 1) * (sizeA + nBasisFuns)
    DO i = 1, nObs
      DO j = 1, nBasisFuns
        As(ithPtStart + i + (j-1) * nObs) = Quadratic(st3d(i), j)
      ENDDO
    ENDDO
    return
  END SUBROUTINE constructiThObsMat

  SUBROUTINE preSolveiThLS(ith, as)
!    use nim_lapack, only: sgeqr2

    INTEGER, intent(in) :: ith
    real, intent(inout) :: as(noa)

    INTEGER :: i, m, n, lda, info, sizeA, ithPtStart
    REAL :: wk(nObs*nBasisFuns)
    real :: work(nBasisFuns)

    m = nObs
    n = nBasisFuns
    lda = m
    sizeA = m * n

    ithPtStart = (ith - 1) * (sizeA + n)
!JR Not sure what is going on with wk(n+1)
    CALL sgeqr2(m, n, As(ithPtStart+1), lda, wk(1), wk(n+1), info)
    DO i = 1, n
      As(ithPtStart + sizeA + i) = wk(i)
    ENDDO

    return
  END SUBROUTINE preSolveiThLS

!---------------------------------------------------------------------------
!  This file contains the classes that implement the 3 variable quadratic 
!  function: 
!
!  Ning Wang, Feb. 2009
!  
!---------------------------------------------------------------------------
  REAL FUNCTION Quadratic(site, bi)
    TYPE(Site3d), intent(in) :: site
    INTEGER, intent(in) :: bi

    SELECT CASE (bi) 
    CASE (1)   
      Quadratic = site%x * site%x

    CASE (2)   
      Quadratic = site%y * site%y
      
    CASE (3)   
      Quadratic = site%z * site%z
      
    CASE (4)   
      Quadratic = site%x * site%y
      
    CASE (5)   
      Quadratic = site%x * site%z
      
    CASE (6)   
      Quadratic = site%y * site%z
      
    CASE (7)   
      Quadratic = site%x
      
    CASE (8)   
      Quadratic = site%y
      
    CASE (9)   
      Quadratic = site%z
      
    CASE (10)   
      Quadratic = 1.
      
    CASE default  
      PRINT*,  'Wrong basis function index! Returning infinity' 
      Quadratic = inf
      
    END SELECT
    
  END FUNCTION Quadratic
end module amtx
