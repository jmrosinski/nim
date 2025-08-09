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

! Compare values in two vectors (avector and bvector) and print difference
! statistics.
! Uses method borrowed from WRF2.2 (see external/io_netcdf/diffwrf.F90).
! JFM added digits2 because digits one was often too large.
subroutine diff_vectors(testName,avector,bvector,npoints,print_header)
  implicit none
  character*(*),intent(in   ) :: testName
  real,         intent(in   ) :: avector(npoints)
  real,         intent(in   ) :: bvector(npoints)
  integer,      intent(in   ) :: npoints
  logical,      intent(in   ) :: print_header
  integer :: i
  ! Difference statistics borrowed from external/io_netcdf/diffwrf.F90
  ! in WRF2.2.
  integer :: ifdiffs, n
  real*8 :: a, b,as,bs,bs2,suma,sumb,meana,meanb
  real*8 :: sumE, sum1, sum2, diff1, diff2, serr, perr, rmse, rms1, rms2, tmp1
  integer digits,d1, d2
  integer digits2,is

  IFDIFFS=0
  sumE  = 0.0
  sum1  = 0.0
  sum2  = 0.0
  suma  = 0.0
  sumb  = 0.0
  diff1 = 0.0
  diff2 = 0.0
  as    = 0.0
  bs    = 0.0
  bs2   = 0.0
  is    = -88
  n = 0
  do i=1,npoints
    a = avector(i)
    b = bvector(i)
    ! borrowed from  Thomas Oppe's comp program
    sumE = sumE + ( a - b ) * ( a - b )
    sum1 = sum1 + a * a
    sum2 = sum2 + b * b
    suma = suma + abs(a)
    sumb = sumb + abs(b)
    if( abs ( a - b ) > diff1) then
      diff1 = abs ( a - b ) 
      is = i
      as = a
      bs = b
    endif
    if( abs ( b ) > diff2) then
      diff2 = abs(b)
      bs2 = b
    endif
    n = n + 1
    IF (a .ne. b) then
      IFDIFFS = IFDIFFS + 1
    ENDIF
  enddo
  meana = suma/npoints
  meanb = sumb/npoints
  rmsE = sqrt ( sumE / dble( n ) )
  rms1 = sqrt ( sum1 / dble( n ) )
  rms2 = sqrt ( sum2 / dble( n ) )
  serr = 0.0
  IF ( sum2 .GT. 0.0d0 ) THEN
    serr = sqrt ( sumE / sum2 )
  ELSE
    IF ( sumE .GT. 0.0d0 ) serr = 1.0
  ENDIF
  perr = 0.0
  IF ( diff2 .GT. 0.0d0 ) THEN
    perr = diff1/diff2
  ELSE
    IF ( diff1 .GT. 0.0d0 ) perr = 1.0
  ENDIF

  IF ( rms1 - rms2 .EQ. 0.0d0 ) THEN
    digits = 8
  ELSE
    IF ( rms2 .NE. 0 ) THEN
      tmp1 = 1.0d0/( ( abs( rms1 - rms2 ) ) / rms2 )
      IF ( tmp1 > 0 ) THEN
        digits = log10(tmp1)
      else
        digits = 9
      ENDIF
    ENDIF
  ENDIF
  if(diff1==0.0) then
    digits2 = 8
  elseif(diff2==0.0) then
    digits2 = 8
  else
    digits2 = min(8,int(abs(log10(diff1/diff2))))
  endif
  IF (print_header) THEN
    PRINT 76
  ENDIF
  digits2=min(digits,digits2)
  IF (IFDIFFS .NE. 0 ) THEN
     PRINT 77,trim(testName), IFDIFFS, rms1, rms2, rmsE, perr,digits2,as,bs,bs2,meana,meanb,is
  ELSE
     PRINT 78,trim(testName), IFDIFFS
  ENDIF
 76 FORMAT (2x,'Variable ',3x,'Ndifs',8x,'RMS (1)',12x,'RMS (2)',6x, &
            'RMSE',5x,'pntwise max',2x,'DIGITS',2x,'pnt1',8x,'pnt2',8x,'pnt2b',7x,'mean1',7x,'mean2')
 77 FORMAT (A10,I9,e18.10,e18.10,e12.4,e12.4,i6,5e12.4,i8 )
 78 FORMAT (A10,I9)
  return
end subroutine diff_vectors
