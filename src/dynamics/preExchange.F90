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

subroutine preExchange(barrierName,exchangeName)
use module_control,only : doBarrier
implicit none
#include <gptl.inc>
character(len=*),intent(IN) :: barrierName,exchangeName
integer                     :: ret

ret = gptlstop  ('MainLoopCompute')
if(doBarrier) then
  ret = gptlstart('PreExchangeBarrier')
  ret = gptlstart(barrierName)
!sms$barrier
  ret = gptlstop (barrierName)
  ret = gptlstop ('PreExchangeBarrier')
endif
ret = gptlstart ('MainLoopExchange')
ret = gptlstart (exchangeName)

return
end subroutine preExchange

