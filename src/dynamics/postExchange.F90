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

subroutine postExchange(barrierName,exchangeName)
use module_control,only : doBarrier
implicit none
#include <gptl.inc>
character(len=*),intent(IN) :: barrierName,exchangeName
integer                     :: ret

ret = gptlstop (exchangeName)
ret = gptlstop ('MainLoopExchange')
if(doBarrier) then
  ret = gptlstart('PostExchangeBarrier')
  ret = gptlstart(barrierName)
!sms$barrier
  ret = gptlstop (barrierName)
  ret = gptlstop ('PostExchangeBarrier')
endif
ret = gptlstart ('MainLoopCompute')

return
end subroutine postExchange

