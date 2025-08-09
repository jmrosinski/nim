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

module module_die
  
contains

  subroutine die(istatus,message,oldstatus)

#ifndef SERIAL
    use mpi
#endif
    
    ! End program with message

    implicit none

    integer,intent(inout)::istatus
    character(len=*),intent(in)::message
    integer,intent(in),optional::oldstatus
    integer::ignore

    if (present(oldstatus)) then
!sms$ignore begin
      write (*,"(a,i0)") trim(message),oldstatus
!sms$ignore end
    else
!sms$ignore begin
      write (*,*) trim(message)
!sms$ignore end
    end if
    call flush(6)
!SMS$INSERT call mpi_abort(mpi_comm_world,istatus,ignore)
    stop

  end subroutine die

  subroutine die_if(condition,istatus,message,oldstatus)

    ! End program and output message if condition is true
    
    implicit none
    logical,intent(in)::condition
    integer,intent(inout)::istatus
    character(len=*),intent(in)::message
    integer,intent(in),optional::oldstatus

    if (condition) then
      if (present(oldstatus)) then
        call die(istatus,message,oldstatus)
      else
        call die(istatus,message)
      end if
    end if

  end subroutine die_if

end module module_die
