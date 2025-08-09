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

subroutine barrier(name,callBarrier)
implicit none
#include <gptl.inc>
character(len=*),intent(IN) :: name
logical         ,intent(IN) :: callBarrier
integer                     :: ret

if(callBarrier) then
             ret = gptlstart('Barrier')
             ret = gptlstart(name)
!sms$barrier
             ret = gptlstop (name)
             ret = gptlstop ('Barrier')
endif

return
end subroutine barrier

