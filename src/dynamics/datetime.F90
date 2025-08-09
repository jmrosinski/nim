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

subroutine datetime
character (24) :: DateT
character ( 8) :: date
character (10) :: time
character ( 5) :: zone
character ( 3) :: month(12)
character (80) :: FMT
integer                    :: values(8)
data month /'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'/

call date_and_time(date,time,zone,values)
FMT = "(a3,i3,',',i5,i3,':',i2.2,':',i2.2)"
write(DateT,FMT) month(values(2)),values(3),values(1),values(5),values(6),values(7)
print "(' DATE-TIME:  ',A)",DateT

end subroutine datetime

!**************************************************************************
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
!**************************************************************************
