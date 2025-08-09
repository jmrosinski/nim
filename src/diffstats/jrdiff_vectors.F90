! Compare values in two vectors (avector and bvector) and print difference
! statistics.

subroutine diff_vectors (testname, avector, bvector, npoints, print_header,CutoffFactor,digitsDmax,BinDigits)
  implicit none

  character(len=*), intent(in) :: testname
  real   , intent(in)  :: avector(npoints)
  real   , intent(in)  :: bvector(npoints)
  integer, intent(in)  :: npoints  ! number of points to compare for each field
  logical, intent(in)  :: print_header
  real   , intent(in)  :: CutoffFactor !Fraction of abs(avector+bvector) below which BinDigits is not counted.
  real   , intent(out) :: digitsDmax !The number of digits to which avector and bvector match
  integer, intent(out) :: BinDigits(0:8) !Histogram of matching digits

  ! Difference statistics borrowed from external/io_netcdf/diffwrf.F90
  ! in WRF2.2.
  real*8 :: a, b               ! point values for case A and B
  real*8 :: rmse               ! RMS error
  real*8 :: delta              ! case A minus case B difference at a point
  real*8 :: deltasave_rdiffmax ! difference for point with biggest relative diff
  real*8 :: rdiff              ! relative difference
  real*8 :: diffmax            ! max pointwise absolute difference
  real*8 :: rdiffmax           ! max pointwise relative difference
  real*8 :: asave_diffmax, bsave_diffmax       ! point values which result in max diff
  real*8 :: asave_rdiffmax, bsave_rdiffmax     ! poing values which result in max relative diff
  real*8 :: digits             ! number of digits that match
  real*8 :: meandigits         ! mean number of digits that match
  real*8 :: digits_diffmax     ! number of digits that match for biggest diff (flag for when no diffs)
  real*8 :: digits_rdiffmax    ! number of digits that match for biggest relative diff (flag no diffs)
  real   :: suma,sumb          ! sum of absolute value of avector and bvector
  real   :: meanA,meanB        ! mean of absolute value of avector and bvector
  real   :: cutoff             ! For calculating BinDigits, only count values above cutoff
  integer:: i,n                ! index
  real*8 :: ndiff              ! number of field values where diffs were observed

  ndiff = 0.
  meandigits = 0.
  digits_diffmax = -999.   ! flag value for no diffs
  digits_rdiffmax = -999.  ! flag value for no diffs
  asave_rdiffmax = -999.
  bsave_rdiffmax = -999.

  diffmax = -999.   ! any negative init value will do
  rdiffmax = -999.  ! any negative init value will do

  suma = 0.0
  sumb = 0.0
  rmse = 0.0
  do i=1,npoints
    suma = suma + abs(avector(i))
    sumb = sumb + abs(bvector(i))
  enddo
  meanA = suma/npoints
  meanB = sumb/npoints
  cutoff = CutoffFactor*(meanA+meanB)/2
  do i=1,npoints
    a = avector(i)
    b = bvector(i)
    delta = abs (a - b)
    rmse = rmse + delta*delta
    if (delta > diffmax) then
      diffmax = delta
      asave_diffmax = avector(i)
      bsave_diffmax = bvector(i)
    end if

    if (delta > 0.) then
      ndiff = ndiff + 1.
      rdiff = delta / (0.5*(abs(a) + abs(b)))
      digits = -log10 (rdiff)                   ! e.g. rdiff=0.1 means digits=1
      if(abs(a+b)/2 > cutoff) then
        n = nint(digits)
        if(n < 0 .or. n > 8) then
          print*,'Error in jrdiff_vectors: nint(digits) out of range',n,digits
        endif
        BinDigits(n) = BinDigits(n) + 1
      endif
      meandigits = meandigits + digits

      if (rdiff > rdiffmax) then
        rdiffmax = rdiff
        asave_rdiffmax = avector(i)
        bsave_rdiffmax = bvector(i)
        deltasave_rdiffmax = delta
      end if
    else
      BinDigits(8) = BinDigits(8) + 1
    end if
  end do

! Mean number of digits
  rmse = sqrt (rmse/npoints)
  meandigits = meandigits/ndiff

  if (ndiff > 0.) then
! Number of digits for worst absolute difference
    rdiff = diffmax / (0.5*(abs(asave_diffmax) + abs(bsave_diffmax)))
    if (rdiff < 1.) then
      digits_diffmax = -log10 (rdiff)  ! e.g. rdiff=0.1 means digits=1
    else
      digits_diffmax = 0
    end if
  
! Number of digits for worst relative difference
    rdiff = deltasave_rdiffmax / (0.5* (abs(asave_rdiffmax) + abs(bsave_rdiffmax)))
    if (rdiff < 1.) then
      digits_rdiffmax = -log10 (rdiff)  ! e.g. rdiff=0.1 means digits=1
    else
      digits_rdiffmax = 0
    end if
  end if

  if (print_header) then
    write(6,*)
    write(6,'(a)')'Notes on difference info printed below:'
    write(6,'(a)')'Max and Min are field maximum and minimum for each field'
    write(6,'(a)')'Mean Abs is the mean of the absolute value of each field'
    write(6,'(a)')'ndiff is the number of points that were not identical between cases'
    write(6,'(a)')'Max diff followed by Values is the max pointwise difference, and field values.'
    write(6,'(a)')'1st digits column is the number of decimal digits that matched for the Max diff point'
    write(6,'(a)')'Max rdiff followed by Values is the max relative difference ((A-B)/avg(A,B)), and field values.'
    write(6,'(a)')'2nd digits column is the number of decimal digits that matched for the Max rdiff point'
    write(6,'(a)')'Avg digits is the average number of digits that matched for all points that were not identical'
    write(6,*)
    print 76
76  format('Field    Max        Min      Mean Abs   ndiff   Max diff      Values     digits  ', &
         'Max rdiff     Values     digits  Avg digits')
  end if

  if (ndiff > 0) then
    write(6,100) testname,maxval(avector),minval(avector),meanA,ndiff,diffmax,asave_diffmax,digits_diffmax, &
                 rdiffmax, asave_rdiffmax, digits_rdiffmax, meandigits
    write(6,101)          maxval(bvector),minval(bvector),meanB,              bsave_diffmax, &
                           bsave_rdiffmax
    digitsDmax = digits_diffmax
  else
    write(6,100) testname, maxval(avector), minval(avector), meanA, ndiff, diffmax
    write(6,101)           maxval(bvector), minval(bvector), meanB
    digitsDmax = 8.
  end if
100 format (a5, 1p, 3e11.3, e8.1, e13.6, e14.6, 0p, f5.1, 1p, e14.6, e14.6, 0p, f5.1, f8.1)
101 format (5x, 1p, 3e11.3, 21x,         e14.6,     19x,             e14.6                )

  return
end subroutine diff_vectors
