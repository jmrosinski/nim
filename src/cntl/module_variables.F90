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
!	This module specifies the variables used in the nim model
!  	Jin-Luen Lee        		December,   2008
!  	Modified from module coded by A. E. MacDonald	October 11, 2004
!	New variable names installed by A. E. MacDonald December 29, 2009
!*********************************************************************

module module_variables
use kinds, only: rt
use module_constants
use module_control
implicit none
save

!.....................................................................
!	Sec. 0.  Dimension Parameters
!.....................................................................
integer	, parameter :: nx=64	!  rhombus x dimension
integer	, parameter :: ny=64	!  rhombus y dimension
integer	, parameter :: nr=10	!  "number of rhombi"
integer	, parameter :: npp=6	!  number of proximity points(max)

!.....................................................................
!	Sec. 1.  Index Variables
!.....................................................................
!  Loop indices
integer its		!  index time step

! Adams Bashforth field indices
integer irk     ! index RK step

!.....................................................................
!	Sec. 2.  Grid Specification Arraies
!.....................................................................

real(rt),allocatable :: etac(:)   ! terrain following mass coord, eta-layer
real(rt),allocatable :: etae(:)   ! terrain following mass coord, eta-edge
real(rt),allocatable :: etin(:)   ! etae for hydrostatic mean data
real(rt),allocatable :: pm1d(:)   ! mean pressure
real(rt),allocatable :: tm1d(:)   ! mean temperature
real(rt),allocatable :: rm1d(:)   ! mean density
real(rt),allocatable :: th1d(:)   ! mean potential temperature
real(rt),allocatable :: zh1d(:)   ! mean height
real(rt) :: pmax,pmin             ! max and min scalar for fct


!....................................................................
!	Sec. 13. Logical variables
!....................................................................
!integer full_ts                  ! logical variable for "full time step" (1=yes)
logical             :: flgrestart ! logical variable for restart

!....................................................................
!	Sec. 14.  Other Variables
!....................................................................
real(rt)   ,allocatable :: ca4k    (:)
real(rt)   ,allocatable :: ca4p    (:)
real(rt)   ,allocatable :: caw4tk  (:)
real(rt)   ,allocatable :: caw4tp  (:)
real(rt)   ,allocatable :: ct4w    (:)

!  dissipation coefs
real(rt) cfh                          ! horizontal numerical dissipation coef
integer im1,ip1		          ! i minus 1, i plus 1
real(rt),allocatable :: spncef(:)     ! spnge layers coefficient
real(rt) ,allocatable :: diagt1d(:,:)      ! time series diagnostic variables

real(rt), ALLOCATABLE :: A(:)

integer,allocatable :: globalperm(:) ! global array containing the permutation of the grid
!NOTE:  TARGET attribute required by icosio!!  
integer,allocatable,target :: global_inv_perm(:) ! inverse permutation of the grid

!SMS$DISTRIBUTE(dh,1) BEGIN
integer,allocatable :: nprox   (:) ! Holds number of proximity points
integer,allocatable ::     perm(:) ! permutation of the grid
real(rt),allocatable :: zs      (:)
!real(rt),allocatable :: r_plus  (:) ! Zaleseks r plus
!real(rt),allocatable :: r_mnus  (:) ! Zaleseks r minus
real(rt),allocatable :: area    (:) ! The area of cell polygon
real(rt),allocatable :: f       (:) ! Coriolis acceleration
real(rt),allocatable :: lat     (:) ! Latitude in radians
real(rt),allocatable :: lon     (:) ! Longitude in radians
real(rt),allocatable :: deg_lat (:) ! lat in degrees
real(rt),allocatable :: deg_lon (:) ! lon in degrees
real(rt),allocatable :: map2    (:) ! (map factor)**2
real(rt),allocatable :: zbs     (:)
real(rt),allocatable :: zbsx    (:)
real(rt),allocatable :: zb      (:)
real(rt),allocatable :: zbx     (:)
real(rt),allocatable :: st2d    (:) ! sea surface temperature

! Diagnostic variables for physics
real(rt),allocatable :: rc2d        (:)  ! convective-scale precipitation
real(rt),allocatable :: rn2d        (:)  ! nonconvective-scale precipitation
real(rt),allocatable :: sn2d        (:)  ! snow from WSM3
real(rt),allocatable :: sw2d(:)          ! downward short-wave radiation flux
real(rt),allocatable :: lw2d(:)          ! downward long-wave radiation flux
real(rt),allocatable :: sfcr_dswF   (:)  ! SFC downward SW flux
real(rt),allocatable :: sfcr_dlwF   (:)  ! SFC downward LW flux
real(rt),allocatable :: sfc_shflx   (:)  ! SFC sensible heat flux
real(rt),allocatable :: sfc_lhflx   (:)  ! SFC latent heat flux
real(rt),allocatable :: toarF       (:)  ! TOA downward SW flux
real(rt),allocatable :: olrF        (:)  ! Outgoint LW radiation flux
!SMS$DISTRIBUTE END

!SMS$DISTRIBUTE(dh,2) BEGIN
! Proximity grid indices:
real(rt),allocatable :: zbxb (:,:)
integer,allocatable  :: prox (:,:)  !  holds index of proximity points
integer,allocatable  :: proxs(:,:)  !  holds index of proximity sides
real(rt),allocatable :: zsc  (:,:)
real(rt),allocatable :: zsb  (:,:)
real(rt),allocatable :: z    (:,:)
real(rt),allocatable :: zm   (:,:)
real(rt),allocatable :: u    (:,:)  ! zonal wind (m/s)
real(rt),allocatable :: v    (:,:)  ! meridional wind (m/s)
real(rt),allocatable :: w    (:,:)  ! vertical velocity (m/s)
real(rt),allocatable :: r    (:,:)  ! density (kg/m**3) r for "rho"
real(rt),allocatable :: t    (:,:)  ! potential temperature (K)
real(rt),allocatable :: tkv  (:,:)  ! virtual temperature (K)
real(rt),allocatable :: ur   (:,:)  ! u*r, flux variable u times density
real(rt),allocatable :: vr   (:,:)  ! v*r, flux variable v times density
real(rt),allocatable :: wr   (:,:)  ! w*r, flux variable w times density
real(rt),allocatable :: tr   (:,:)  ! t*r, flux variable t times density
real(rt),allocatable :: p    (:,:)  ! pressure (pascals)
real(rt),allocatable :: e    (:,:)  ! Exner func = (p/p0)**kappa
real(rt),allocatable :: eb   (:,:)  ! Exner func = (p/p0)**kappa
real(rt),allocatable :: rb   (:,:)  ! density basic state (kg/m**3)
real(rt),allocatable :: tb   (:,:)  ! potential temperature basic state (k)
real(rt),allocatable :: pb   (:,:)  ! pressure basic state, hydrostatic(pascals)
real(rt),allocatable :: trb  (:,:)  ! flux potential temperature basic state
real(rt),allocatable :: rp   (:,:)  ! rho (density) prime (kg/m**3) (= r - rb)
real(rt),allocatable :: tp   (:,:)  ! potential temperature prime (k) (= t - tb)
real(rt),allocatable :: trp  (:,:)  ! flux pot temp prime ( = tr - trb)	
real(rt),allocatable :: qv   (:,:)  ! water vapor specific humidity (kg/kg)
real(rt),allocatable :: qw   (:,:)  ! cloud water (kg/kg of moist air)
!real(rt),allocatable :: qi   (:,:)  ! cloud ice (kg/kg)
real(rt),allocatable :: qr   (:,:)  ! rain (kg/kg)
!real(rt),allocatable :: qs   (:,:)  ! snow (kg/kg)
!real(rt),allocatable :: qg   (:,:)  ! graupel (kg/kg) 
real(rt),allocatable :: rh   (:,:)  ! relative humidity from 0 to 1.
real(rt),allocatable :: div  (:,:)  ! 3D divergence out of cell
real(rt),allocatable :: cs   (:,:)  ! cosine transform function
real(rt),allocatable :: sn   (:,:)  ! sine transform function
real(rt),allocatable :: map  (:,:)  ! map factor at mid point of edge chords
real(rt),allocatable :: rmap (:,:)  ! reciprocal of "map(6,nip)"
real(rt),allocatable :: zbse (:,:)  ! topo on the mid-pt of sidevec
real(rt),allocatable :: zbe  (:,:)  ! topo on the mid-pt of sidevec
real(rt),allocatable :: zclvl(:,:)  ! z at level on the center of the cell
real(rt),allocatable :: vol  (:,:)  ! volume for each control volume
real(rt),allocatable :: su   (:,:)  ! source-sink for u equation
real(rt),allocatable :: sv   (:,:)  ! source-sink for v equation
real(rt),allocatable :: sw   (:,:)  ! source-sink for w equation
real(rt),allocatable :: st   (:,:)  ! source-sink for theta equation
real(rt),allocatable :: he3d (:,:)  ! 3-D Gauss-function for heating function
real(rt),allocatable :: tfur (:,:)  ! tendency function of ur
real(rt),allocatable :: tfvr (:,:)  ! tendency function of vr
real(rt),allocatable :: tfwr (:,:)  ! tendency function of wr
real(rt),allocatable :: tftr (:,:)  ! tendency function of tr
real(rt),allocatable :: tfr  (:,:)  ! tendency function of density (r)
real(rt),allocatable :: fur  (:,:)  ! forcing of ur
real(rt),allocatable :: fvr  (:,:)  ! forcing of vr
real(rt),allocatable :: fwr  (:,:)  ! forcing of wr
real(rt),allocatable :: ftr  (:,:)  ! forcing of tr
real(rt),allocatable :: fr   (:,:)  ! forcing of density (r)
real(rt),allocatable :: zb3d (:,:)
real(rt),allocatable :: rebb (:,:)  ! rb at edge
real(rt),allocatable :: trebb(:,:)  ! trb at edge
real(rt),allocatable :: tebb (:,:)  ! tb at edge

real(rt),allocatable  :: fu0  (:,:)
real(rt),allocatable  :: fv0  (:,:)
real(rt),allocatable  :: fw0  (:,:)

real(rt),allocatable ::  rs  (:,:)  
real(rt),allocatable :: urs  (:,:)  
real(rt),allocatable :: vrs  (:,:)  
real(rt),allocatable :: wrs  (:,:)  
real(rt),allocatable :: trs  (:,:)  
real(rt),allocatable :: fdt  (:)    ! fraction of rk4 time step
real(rt),allocatable :: frk  (:)    ! factor of rk4 increments

real(rt), allocatable :: nvecb  (:,:,:)
real(rt), allocatable :: vdnb   (:,:,:)
real(rt), allocatable :: uw8b   (:,:,:) !Bu(0:nz),Bv(0:nz),Bw(0:nz)
real(rt), allocatable :: bedgvar(:,:,:)
real(rt),allocatable :: trq1  (:,:,:)  ! tracers group 1 : trq1=trc1/r
real(rt),allocatable :: sedgtrc(:,:,:,:)
real(rt),allocatable :: bedgtrc(:,:,:)

! Diagnostic variables for physics
real(rt),allocatable :: atmr_lwH    (:,:)  ! ATM LW heating rate
real(rt),allocatable :: atmr_swH    (:,:)  ! ATM SW heating rate
real(rt),allocatable :: atm_turbH   (:,:)  ! ATM turbulent heat flux

!SMS$DISTRIBUTE END

!SMS$DISTRIBUTE(dh,3) BEGIN
real(rt),allocatable :: trc1   (:,:,:)  ! tracers group 1 : 1:r*qv, 2:r*qw, 3:O3
real(rt),allocatable :: xyc    (:,:,:)
real(rt),allocatable :: xys    (:,:,:)
real(rt),allocatable :: xyb    (:,:,:)
real(rt),allocatable :: sa     (:,:,:)
real(rt),allocatable :: zc     (:,:,:)
real(rt),allocatable :: zcs    (:,:,:)
real(rt),allocatable :: reb    (:,:,:)  ! rb at edge
real(rt),allocatable :: ue     (:,:,:)  ! u at edge, center point of chord
real(rt),allocatable :: ve     (:,:,:)  ! v at edge, center point of chord
real(rt),allocatable :: we     (:,:,:)  ! w at edge, center point of chord
real(rt),allocatable :: re     (:,:,:)  ! density at edge, center point of chord
real(rt),allocatable :: treb   (:,:,:)  ! trb at edge
real(rt),allocatable ::  teb   (:,:,:)  ! tb at edge
real(rt),allocatable :: ueb    (:,:,:)  ! u at edge, center point of chord
real(rt),allocatable :: veb    (:,:,:)  ! v at edge, center point of chord
real(rt),allocatable :: ure    (:,:,:)  ! ur at edge, center point of chord
real(rt),allocatable :: vre    (:,:,:)  ! vr at edge, center point of chord
real(rt),allocatable :: wre    (:,:,:)  ! wr at edge, center point of chord
real(rt),allocatable :: tre    (:,:,:)  ! tr at edge, center point of chord
real(rt),allocatable :: trpe   (:,:,:)  ! trp at edge, center point of chord
real(rt),allocatable :: rebs   (:,:,:)
real(rt),allocatable :: urebs  (:,:,:)
real(rt),allocatable :: trebs  (:,:,:)
real(rt),allocatable :: trpebs (:,:,:)
real(rt),allocatable :: wrebs  (:,:,:)
real(rt),allocatable :: rep    (:,:,:)  ! rho (density) at edge center, prime
real(rt),allocatable :: tefr   (:,:,:)  ! edge tendency function of density (r)
real(rt),allocatable :: efr    (:,:,:)  ! edge forcing of density (r)
real(rt),allocatable :: prox_xy(:,:,:)  ! holds x and y locations for prox pts
real(rt),allocatable :: zelvl  (:,:,:)  ! z at level on the mid-pt of sidevec
real(rt),allocatable :: saoprx (:,:,:)  
real(rt),allocatable :: vdns   (:,:,:)
real(rt),allocatable :: nvecs  (:,:,:,:)
real(rt),allocatable :: amtx1  (:,:,:)  ! layer variables
real(rt),allocatable :: amtx2  (:,:,:)  ! level variables
real(rt),allocatable :: stmp   (:,:,:)
real(rt),allocatable :: btmp   (:,:,:)

!real(rt),allocatable :: hfxn   (:,:,:)
!real(rt),allocatable :: hfxp   (:,:,:)
!
!SMS$DISTRIBUTE END

!SMS$DISTRIBUTE(dh,4) BEGIN
real(rt),allocatable :: uw8s   (:,:,:,:)!Su(1:nz),Sv(1:nz)
real(rt),allocatable :: sedgvar(:,:,:,:)
!SMS$DISTRIBUTE END

end module module_variables

