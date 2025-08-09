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

module module_alloc
!*********************************************************************
!
!  	Jacques Middlecoff     		March 2009
!
!*********************************************************************
contains
subroutine alloc

use infnan, only: inf
use module_constants
use module_control
use module_variables
use ReadNamelist, only: nz
!use print_taskinfo, only: printmem

use kinds, only: rt
implicit none
#include <gptl.inc>

save
integer :: ipn
integer :: resultlen,ierr,ret
    integer :: procsiz, rss, share, text, datastack     ! returned by gptlget_memusage                                                                                                                                                       

!.....................................................................
!	Sec. 2.  Grid Specification Arrays
!.....................................................................
!

allocate(zs(ims:ime),zsc(npp,ims:ime),zsb(npp,ims:ime))
allocate(xyc(2,npp,ims:ime))
allocate(xys(2,npp,ims:ime))
allocate(xyb(2,npp,ims:ime))
allocate(z(0:nz,ims:ime))
allocate(zm(nz,ims:ime))
allocate(zc(0:nz,npp,ims:ime))
allocate(zcs(nz  ,npp,ims:ime))
allocate(zb3d(0:nz,ims:ime))
allocate(nvecs(nz,npp,ims:ime,2),nvecb(0:nz,ims:ime,3))
allocate(sa(nz,npp,ims:ime))
allocate(etac(1:nz))   ! terrain following mass coord, eta-layer
allocate(etae(0:nz))   ! terrain following mass coord, eta-edge
allocate(etin(1:nz))   ! etae for hydrostatic mean data
allocate(pm1d(0:nz))   ! mean pressure
allocate(tm1d(0:nz))   ! mean temperature
allocate(rm1d(0:nz))   ! mean density
allocate(th1d(0:nz))   ! mean potential temperature
allocate(zh1d(0:nz))   ! mean height

allocate(st2d(ims:ime))

! Su(0) used in W-eq, Bur(0) used in U-eq for vertical adv
allocate(uw8s(1:nz,npp,2,ims:ime)) !Su(1:nz),Sv(1:nz)
allocate(uw8b(0:nz,ims:ime,3)) !Bu(0:nz),Bv(0:nz),Bw(0:nz)
allocate(sedgvar(nz,npp,nvars,ims:ime))

if(                  loc(sedgvar) - loc(sedgvar(lbound(sedgvar,1),lbound(sedgvar,2),lbound(sedgvar,3),lbound(sedgvar,4))) /= 0) then
  print*,'JFMalloc' ,loc(sedgvar) , loc(sedgvar(lbound(sedgvar,1),lbound(sedgvar,2),lbound(sedgvar,3),lbound(sedgvar,4))),ims,ime
endif

allocate(bedgvar(0:nz,ims:ime,nvars))

!.....................................................................
!	Sec. 2.  Prognostic State Variables
!.....................................................................
! State variables:
allocate(u(  nz,ims:ime))	! zonal wind (m/s)
allocate(v(  nz,ims:ime))	! meridional wind (m/s)
allocate(w(0:nz,ims:ime))	! vertical velocity (m/s)
allocate(r(  nz,ims:ime))	! density (kg/m**3) r for "rho"
allocate(t(  nz,ims:ime))	! potential temperature (K)
allocate(tkv(nz,ims:ime))	! air virtual temperature (K)

! Flux variables:
allocate(ur(  nz,ims:ime))	! u*r, flux variable u times density
allocate(vr(  nz,ims:ime))	! v*r, flux variable v times density
allocate(wr(0:nz,ims:ime))	! w*r, flux variable w times density
allocate(tr(  nz,ims:ime))	! t*r, flux variable t times density

!allocate(hfxn(nz,npp,ims:ime))
!allocate(hfxp(nz,npp,ims:ime))

allocate( rs(  nz,ims:ime))	
allocate(urs(  nz,ims:ime))	
allocate(vrs(  nz,ims:ime))	
allocate(wrs(0:nz,ims:ime))	
allocate(trs(  nz,ims:ime))	

allocate(fdt(nrkl)) 
allocate(frk(nrkl)) 

!....................................................................
!	Sec. 3. Diagnostic State Variables
!....................................................................
! Diagnostic main variables:
allocate(p(0:nz,ims:ime))	        ! pressure (pascals)
allocate(e(nz,ims:ime))	        ! Exner func = (p/p0)**kappa
allocate(eb(nz,ims:ime))	        ! Exner func = (p/p0)**kappa

! Basic state variables:
allocate(rb  (  nz    ,ims:ime))	! density basic state (kg/m**3)
allocate(tb  (  nz    ,ims:ime))	! potential temperature basic state (k)
allocate(pb  (0:nz    ,ims:ime))	! pressure basic state, hydrostatic (pascals)
allocate(reb (  nz,npp,ims:ime))	! rb at edge
allocate(rebb(0:nz,ims:ime))	! rb at edge


! Flux variables basic state:
allocate(trb(nz,ims:ime))		! flux potential temperature basic state

! Pertubation main variables:
allocate(rp(nz,ims:ime))		! rho (density) prime (kg/m**3) (= r - rb)
allocate(tp(nz,ims:ime))		! potential temperature prime (k) (= t - tb)

! Perturbation flux variables:
allocate(trp(nz,ims:ime))		! flux pot temp prime ( = tr - trb)	

!......................................................................
!	Sec. 4. Water (microphysical) variables
!......................................................................
allocate(qv(nz,ims:ime))		! water vapor specific humidity (kg/kg)
allocate(qw(nz,ims:ime))		! cloud water (kg/kg of moist air)
!allocate(qi(nz,ims:ime))		! cloud ice (kg/kg)
allocate(qr(nz,ims:ime))		! rain (kg/kg)
!allocate(qs(nz,ims:ime))		! snow (kg/kg)
!allocate(qg(nz,ims:ime))		! graupel (kg/kg) 
allocate(rh(nz,ims:ime))		! relative humidity (0. to 1.)
allocate(trc1(nz,ntr1,ims:ime))	        ! tracers group 1 : 1:r*qv, 2:r*qw, 3:O3
allocate(trq1(nz,ims:ime,ntr1))         ! tracers group 1 : trc1/r
!allocate(sedgtrc(nz,npp,ims:ime,1))    ! tracer values at edges.
!allocate(bedgtrc(nz,ims:ime,0:1))      ! tracer values at edges.


!......................................................................
!	Sec. 5. Edge of cell variables
!......................................................................

! Edge main state variables:
allocate(ue   (  nz,npp,ims:ime))	! u at edge, center point of chord
allocate(ve   (  nz,npp,ims:ime))	! v at edge, center point of chord
allocate(we   (  nz,npp,ims:ime))	! w at edge, center point of chord
allocate(re   (  nz,npp,ims:ime))	! density (rho) at edge, center point of chord
allocate(treb (  nz,npp,ims:ime))   ! trb at edge
allocate(trebb(0:nz,ims:ime))   ! trb at edge
allocate(teb  (  nz,npp,ims:ime))   ! trb at edge
allocate(tebb (0:nz,ims:ime))   ! trb at edge
allocate(ueb  (0:nz,npp,ims:ime))	! u at edge, center point of chord
allocate(veb  (0:nz,npp,ims:ime))	! v at edge, center point of chord

! Edge flux variables:
allocate(ure (  nz,npp,ims:ime))	! ur at edge, center point of chord
allocate(vre (  nz,npp,ims:ime))	! vr at edge, center point of chord
allocate(wre (0:nz,npp,ims:ime))	! wr at edge, center point of chord
allocate(tre (  nz,npp,ims:ime))	! tr at edge, center point of chord
allocate(trpe(  nz,npp,ims:ime))	! trp at edge, center point of chord

allocate( rebs (0:nz,npp,ims:ime))
allocate(urebs (0:nz,npp,ims:ime))
allocate(trebs (0:nz,npp,ims:ime))
allocate(trpebs(0:nz,npp,ims:ime))
allocate(wrebs (0:nz,npp,ims:ime))


! Edge prime variables
allocate(rep(nz,npp,ims:ime))	! rho (density) at edge center, prime

!.....................................................................
!	Sec. 6. Cell calculated variables
!.....................................................................
! Normal velocity used in flux calc in nim:
allocate(vdns(1:nz,npp,ims:ime))
allocate(vdnb(0:nz,ims:ime,2))

! Cell divergence
allocate(div(nz,ims:ime))		! 3D divergence out of cell

! forcing variables:
allocate(tfur(  nz,ims:ime))    ! tendency function of ur
allocate(tfvr(  nz,ims:ime))    ! tendency function of vr
allocate(tfwr(0:nz,ims:ime))    ! tendency function of wr
allocate(tftr(  nz,ims:ime))    ! tendency function of tr
allocate(tfr (  nz,ims:ime))    ! tendency function of density (r)
allocate(tefr (nz,npp,ims:ime)) ! edge tendency function of density (r)
allocate( fur(  nz,ims:ime))    ! forcing of ur
allocate( fvr(  nz,ims:ime))    ! forcing of vr
allocate( fwr(0:nz,ims:ime))    ! forcing of wr
allocate( ftr(  nz,ims:ime))    ! forcing of tr
allocate( fr (  nz,ims:ime))    ! forcing of density (r)
allocate( efr (nz,npp,ims:ime)) ! edge forcing of density (r)
allocate(fu0(nz,ims:ime),fv0(nz,ims:ime),fw0(nz,ims:ime))

! flux conserving transport variables:
!allocate(r_plus(ims:ime),r_mnus(ims:ime))    ! Zaleseks r plus and r minus


!.....................................................................
!	Sec. 8.  Proximity Grid Description Variables
!.....................................................................
! Coordinate tranformation constants - from ll to xy
allocate(cs(npp,ims:ime),sn(npp,ims:ime)) ! cosine and sine transform functions

! Variables to describe the icos grid in xy space
allocate(prox_xy(npp,2,ims:ime))
allocate(saoprx(nz,npp,ims:ime))
allocate(area(ims:ime))	        ! the area of cell polygon

!...................................................................
!	Sec. 9. Van Der Monde interpolation intermediate variables
!...................................................................
! Vandermonde main arrays:

!call printmem ('alloc before amtx1')
!call flush(6)
allocate(amtx1 (  nz  ,(nob+1)*nbf,ims:ime)) ! layer variables
allocate(amtx2 (0:nz-1,(nob+1)*nbf,ims:ime)) ! level variables

!.....................................................................
!	Sec. 10.  Geo Variables:
!.....................................................................
allocate(f      (  ims:ime))  ! Coriolis acceleration
allocate(lat    (  ims:ime))  ! latitude  in radians
allocate(lon    (  ims:ime))  ! longitude in radians
allocate(deg_lat(  ims:ime))  ! latitude  in degrees
allocate(deg_lon(  ims:ime))  ! longitude in degrees
allocate(map    (6,ims:ime))  ! map factor at mid point of edge chords
allocate(rmap   (6,ims:ime))  ! reciprocal of "map(6,ims:ime)"
allocate(map2   (  ims:ime))  ! (map factor)**2

!....................................................................
!       Sec. 11.  Topography
!....................................................................
allocate(zbs(ims:ime),zbsx(ims:ime),zb(ims:ime),zbx(ims:ime),zbxb(npp,ims:ime))

allocate( zbse (     npp,ims:ime))	! topo on the mid-pt of sidevec
allocate( zbe  (     npp,ims:ime))	! topo on the mid-pt of sidevec
allocate( zclvl(0:nz,    ims:ime))	! z at level on the center of the cell
allocate( zelvl(0:nz,npp,ims:ime))	! z at level on the mid-pt of sidevec
allocate( vol  (  nz,    ims:ime))	! volume for each control volume

allocate(ca4k(1:nz))
allocate(ca4p(1:nz))
allocate(caw4tk(1:nz))
allocate(caw4tp(1:nz))
allocate(ct4w(1:nz))

!....................................................................
!	Sec. 12. Source-sink functions
!....................................................................
allocate(su(nz,ims:ime))		! source-sink for u equation
allocate(sv(nz,ims:ime))		! source-sink for v equation
allocate(sw(nz,ims:ime))		! source-sink for w equation
allocate(st(nz,ims:ime))		! source-sink for theta equation

! Set to Infinity all allocated floating point variables. 
! Dont do this in the Lahey case because that compiler has its own way of checking
! referenced but uninitialized variables.


!....................................................................
!	Sec. 13. Physics diagnostic variables 
!....................................................................
allocate( rn2d(ims:ime)      )      ! total precipitation
allocate( rc2d(ims:ime)      )      ! convective precipitation
allocate( sn2d(ims:ime)      )      ! snow from WSM3
allocate( sw2d(ims:ime)      )      ! surface short-wave radiational flux
allocate( lw2d(ims:ime)      )      ! surface long-wave radiational flux
allocate( sfcr_dswF(ims:ime) )      ! SFC downward SW flux
allocate( sfcr_dlwF(ims:ime) )      ! SFC downward LW flux
allocate( sfc_shflx(ims:ime) )      ! SFC sensible heat flux
allocate( sfc_lhflx(ims:ime) )      ! SFC latent heat flux
allocate( toarF    (ims:ime) )      ! TOA downward SW flux
allocate( olrF     (ims:ime) )      ! Outgoint LW radiation flux
allocate( atmr_lwH (nz,ims:ime) )   ! ATM LW heating rate
allocate( atmr_swH (nz,ims:ime) )   ! ATM SW heating rate
allocate( atm_turbH(nz,ims:ime) )   ! ATM turbulent heat flux

!....................................................................
!	Sec. 14.  Other Variables
!....................................................................
allocate(spncef(0:nz))          ! spnge layers coefficient

!TODO:  make this optional via a namelist flag
!allocate(diagt1d(0:100000,10))     ! time series diagnostic variables

! Floating-point variables allocated above should all be set to "inf" 
! below to trigger floating-point exceptions if uninitialized values 
! are used (and compiler flags are set to halt on exceptions).  
! Initialization of variables to other values should be handled 
! in init().  

! Note that we do not set variables to "inf" when using the Lahey compiler 
! because it supports a more sophisticated mechanism for detecting 
! uninitialized values.  

#ifdef LAHEY
rn2d      = 0._rt
rc2d      = 0._rt

#else

!JR Use loops to initialize for optimal first-touch in OMP mode

! this array is not decomposed
!diagt1d(:,:) = inf

!$OMP PARALLEL DO
do ipn=ims,ime
  zs(ipn) = inf
  zsc(:,ipn) = inf
  zsb(:,ipn) = inf
  xyc(:,:,ipn) = inf
  xys(:,:,ipn) = inf
  xyb(:,:,ipn) = inf
  z(:,ipn) = inf
  zm(:,ipn) = inf
  zc(:,:,ipn) = inf
  zcs(:,:,ipn) = inf
  zb3d(:,ipn) = inf
  nvecs(:,:,ipn,:) = inf
  nvecb(:,ipn,:) = inf
  sa(:,:,ipn) = inf
end do

etac(:) = inf
etae(:) = inf
etin(:) = inf
pm1d(:) = inf
tm1d(:) = inf
rm1d(:) = inf
th1d(:) = inf
zh1d(:) = inf

!$OMP PARALLEL DO
do ipn=ims,ime
  st2d(ipn) = inf
  uw8s(:,:,:,ipn) = inf
  uw8b(:,ipn,:) = inf
  sedgvar(:,:,:,ipn) = inf
  bedgvar(:,ipn,:) = inf
  u(:,ipn) = inf
  v(:,ipn) = inf
  w(:,ipn) = inf
  r(:,ipn) = inf
  t(:,ipn) = inf
  tkv(:,ipn) = inf
  ur(:,ipn) = inf
  vr(:,ipn) = inf
  wr(:,ipn) = inf
  tr(:,ipn) = inf
  rs(:,ipn) = inf
  urs(:,ipn) = inf
  vrs(:,ipn) = inf
  wrs(:,ipn) = inf
  trs(:,ipn) = inf
end do

fdt(:) = inf
frk(:) = inf

!$OMP PARALLEL DO
do ipn=ims,ime
  p(:,ipn) = inf
  e(:,ipn) = inf
  eb(:,ipn) = inf
  rb(:,ipn) = inf
  tb(:,ipn) = inf
  pb(:,ipn) = inf
  reb(:,:,ipn) = inf
  rebb(:,ipn) = inf
  trb(:,ipn) = inf
  rp(:,ipn) = inf
  tp(:,ipn) = inf
  trp(:,ipn) = inf
  qv(:,ipn) = inf
  qw(:,ipn) = inf
  rh(:,ipn) = inf
  trc1(:,:,ipn) = inf
  trq1(:,ipn,:) = inf
  ue(:,:,ipn) = inf
  ve(:,:,ipn) = inf
  we(:,:,ipn) = inf
  re(:,:,ipn) = inf
  treb(:,:,ipn) = inf
  trebb(:,ipn) = inf
  teb(:,:,ipn) = inf
  tebb(:,ipn) = inf
  ueb(:,:,ipn) = inf
  veb(:,:,ipn) = inf
  ure(:,:,ipn) = inf
  vre(:,:,ipn) = inf
  wre(:,:,ipn) = inf
  tre(:,:,ipn) = inf
  trpe(:,:,ipn) = inf
  rebs(:,:,ipn) = inf
  urebs(:,:,ipn) = inf
  trebs(:,:,ipn) = inf
  trpebs(:,:,ipn) = inf
  wrebs(:,:,ipn) = inf
  rep(:,:,ipn) = inf
  vdns(:,:,ipn) = inf
  vdnb(:,ipn,:) = inf
  div(:,ipn) = inf
  tfur(:,ipn) = inf
  tfvr(:,ipn) = inf
  tfwr(:,ipn) = inf
  tftr(:,ipn) = inf
  tfr(:,ipn) = inf
  tefr(:,:,ipn) = inf
  fur(:,ipn) = inf
  fvr(:,ipn) = inf
  fwr(:,ipn) = inf
  ftr(:,ipn) = inf
  fr(:,ipn) = inf
  efr(:,:,ipn) = inf
!  r_plus(ipn) = inf
!  r_mnus(ipn) = inf
  cs(:,ipn) = inf
  sn(:,ipn) = inf
  prox_xy(:,:,ipn) = inf
  saoprx(:,:,ipn) = inf
  area(ipn) = inf
  amtx1(:,:,ipn) = inf
  amtx2(:,:,ipn) = inf
  f(ipn) = inf
  lat(ipn) = inf
  lon(ipn) = inf
  deg_lat(ipn) = inf
  deg_lon(ipn) = inf
  map(:,ipn) = inf
  rmap(:,ipn) = inf
  map2(ipn) = inf
  zbs(ipn) = inf
  zbsx(ipn) = inf
  zb(ipn) = inf
  zbx(ipn) = inf
  zbxb(:,ipn) = inf
  zbse(:,ipn) = inf
  zbe(:,ipn) = inf
  zclvl(:,ipn) = inf
  zelvl(:,:,ipn) = inf
  vol(:,ipn) = inf
end do

ca4k(:) = inf
ca4p(:) = inf
caw4tk(:) = inf
caw4tp(:) = inf
ct4w(:) = inf
spncef(:) = inf

!$OMP PARALLEL DO
do ipn=ips,ipe
  su(:,ipn) = inf
  sv(:,ipn) = inf
  sw(:,ipn) = inf
  st(:,ipn) = inf
  rn2d(ipn)      = inf
  rc2d(ipn)      = inf
  sfcr_dswF(ipn) = inf
  sfcr_dlwF(ipn) = inf
  sfc_shflx(ipn) = inf
  sfc_lhflx(ipn) = inf
  sn2d(ipn)      = inf
  sw2d(ipn)      = inf
  lw2d(ipn)      = inf
  toarF(ipn)     = inf
  olrF(ipn)      = inf
  atmr_lwH(:,ipn)  = inf
  atmr_swH(:,ipn)  = inf
  atm_turbH(:,ipn) = inf
end do
#endif

end subroutine alloc

end module module_alloc

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
