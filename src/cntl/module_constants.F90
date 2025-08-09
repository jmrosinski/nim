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

!********************************************************************
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
!       This module specifies model constants 
!       A. E. MacDonald         October 11, 2004
!       J. LEE                  September,  2005
!********************************************************************

MODULE module_constants
use kinds, only: rt
implicit none

save
real(rt), parameter :: pi       = 3.14159263_rt
real(rt), PARAMETER :: DEGRAD   = pi/180._rt
real(rt), parameter :: raddeg   = 180._rt/pi
real(rt), parameter :: g        = 9.80616_rt
real(rt), parameter :: ae       = 6371220._rt   ! earth radius
real(rt), parameter :: omega    = .00007292_rt  ! earth rotation rate

real(rt), parameter :: p1000    = 100000._rt    ! p at 1000mb (pascals)
real(rt), parameter :: cp       = 1004._rt      ! specific heat at const pres
real(rt), parameter :: cv       =  717._rt      ! specific heat at const vol  
real(rt), parameter :: rd       = 287.0586_rt   ! spec gas constant for dry air
real(rt), parameter :: rv       = 461.50_rt     ! gas constant for H2O
real(rt), parameter :: gamma     = cp/cv
real(rt), parameter :: kapa     = rd/cp
real(rt), parameter :: kapi     = rd/cv         ! exponential for Exnr
real(rt), parameter :: dcef     = 0._rt

integer,parameter :: double_real = 0, direct_norm = 0
 
END MODULE module_constants
