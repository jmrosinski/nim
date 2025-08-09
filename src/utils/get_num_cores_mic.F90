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

program get_num_cores_mic
  use ReadNamelist, only: readnl, mpipn_mic, nthreads_mic
  implicit none

  integer :: ret

  call readnl (ret)
  if (ret /= 0) then
    write(6,*)'get_num_cores_mic: readnl gave a bad return code=', ret
    stop 999
  end if

! Check for valid input

  if (mpipn_mic < 1) then
    write (6,*)'get_num_cores_mic: mpipn_mic=', mpipn_mic, ' invalid or not found'
    stop 999
  end if

  if (nthreads_mic < 1) then
    write (6,*)'get_num_cores_mic: nthreads_mic=', nthreads_mic, ' invalid or not found'
    stop 999
  end if

  write (*,'(a,i0)') 'mpipn_mic:', mpipn_mic
  write (*,'(a,i0)') 'nthreads_mic:', nthreads_mic
end program get_num_cores_mic
