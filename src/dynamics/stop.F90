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

subroutine stop
!*********************************************************************
!
!	Stop program for NIM
!	Alexander E. MacDonald  12/27/2004
!
!*********************************************************************

use print_taskinfo, only: printmem
use module_control, only: preExchangeBarrier,postExchangeBarrier
use core_setup    , only: iam_nim_task
use globals       , only: myrank, npes

implicit none

#include <gptl.inc>

integer :: ret,comm
real(8) :: TOTALTIME,MLT,MLE,MLC,EXB ! returned from gptl
real(8) :: tmin, tmax
!SMS$SERIAL BEGIN
close(29)
!SMS$SERIAL END

if (iam_nim_task) then

!sms$unstructured_print_timers

  print*,''
  print*,'                            ROUTINE   MIN TIME(sec)  MAX TIME(sec)'
  ret = gptlget_wallclock ('MainLoop', 0, MLT )
  tmin = MLT
  tmax = MLT
!SMS$REDUCE(tmin,min)
!SMS$REDUCE(tmax,max)
  print"(1x,A35,F13.3,F15.3)", 'Main Loop', tmin, tmax

  if(preExchangeBarrier) then
    ret = gptlget_wallclock ('PreExchangeBarrier' , 0, EXB)
    tmin = EXB
    tmax = EXB
!SMS$REDUCE(tmin,min)
!SMS$REDUCE(tmax,max)
    print"(1x,A35,F13.3,F15.3)", 'Main Loop Pre-Exchange Barrier', tmin,tmax
  endif

  ret = gptlget_wallclock ('MainLoopCompute', 0, MLC)
  tmin = MLC
  tmax = MLC
!SMS$REDUCE(tmin,min)
!SMS$REDUCE(tmax,max)
  print"(1x,A35,F13.3,F15.3)", 'Main Loop Compute', tmin, tmax
  ret = gptlget_wallclock ('MainLoopExchange', 0, MLE)
  tmin = MLE
  tmax = MLE
!SMS$REDUCE(tmin,min)
!SMS$REDUCE(tmax,max)
  print"(1x,A35,F13.3,F15.3)", 'Main Loop Exchange', tmin,tmax

  if(postExchangeBarrier) then
    ret = gptlget_wallclock ('PostExchangeBarrier' , 0, EXB)
    tmin = EXB
    tmax = EXB
!SMS$REDUCE(tmin,min)
!SMS$REDUCE(tmax,max)
    print"(1x,A35,F13.3,F15.3)", 'Main Loop Post-Exchange Barrier', tmin,tmax
  endif

  ret = gptlget_wallclock ('MainLoop_nooutput', 0, MLT)
  tmin = MLT
  tmax = MLT
!SMS$REDUCE(tmin,min)
!SMS$REDUCE(tmax,max)
  print"(1x,A35,F13.3,F15.3)", 'MainLoop_nooutput', tmin, tmax
  print*,' '

  call printmem ('End of run')
  ret = gptlstop ('Total')
  ret = gptlget_wallclock ('Total', 0, TOTALTIME)  ! The "0" is thread number
  ! write detailed "timing.$myrank" file only from root for now
  ! remove this "if" statement to write "timing.$myrank" from every MPI task
  if (myrank == 0 .or. myrank == npes-1) then
    ret = gptlpr (myrank)
  endif
!sms$get_communicator(comm)
!sms$insert ret = gptlpr_summary (comm)
  print*,''
  print*,'Total time =' , TOTALTIME

endif

return
end subroutine stop

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
