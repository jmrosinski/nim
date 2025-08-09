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

module print_taskinfo
#ifndef SERIAL
  use mpi
#endif

  implicit none
  public

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! print_nodeinfo: Print info about rank and node name on which we're running
!                 In serial mode the node name is just "serial"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine print_nodeinfo (iam_nim_task,iam_write_task)
    use globals, only: myrank, mynode, nPEs
    implicit none
    logical, intent(in) :: iam_nim_task,iam_write_task

    integer :: resultlen   ! length return from MPI routines
    integer :: ierr=0     ! error return from MPI routines
    integer :: ios        ! error return from open
    integer :: unitno=31  ! unit number for the mpirank.txt file
    integer :: rank       ! the rank of the MPI process
!TODO:  MPI communicators are not guaranteed to be INTEGER, fix this
    integer :: MPI_Comm   ! the MPI communicator from SMS
    character(len=80) :: tasktype
!SMS$INSERT character(len=MPI_MAX_PROCESSOR_NAME) :: nodeName(0:nPEs-1)

    if (iam_nim_task) then
      tasktype='compute'
    else if (iam_write_task) then
      tasktype='write'
    else
      tasktype='donothing'
    endif

!sms$get_communicator(mpi_comm)
!SMS$INSERT call mpi_gather(mynode  ,MPI_MAX_PROCESSOR_NAME,MPI_character, &
!SMS$INSERT nodeName,MPI_MAX_PROCESSOR_NAME,MPI_character,0,MPI_Comm,ierr)

!SMS$SERIAL (default=ignore) BEGIN
open (unitno, file='mpirank'//trim(tasktype)//'.txt', form="formatted", action='write', status='new', iostat=ios)
if (ios /= 0) then
  write(6,*) 'taskinfo: failure to open output file mpirank.txt - write to stdout'
  unitno=6
end if
#ifdef SERIAL
  write(unitno,*)'Serial run: node=', trim(mynode)
#else
  do rank=0,nPEs-1
    write(unitno,100)'MPI rank=', rank,' is a ',trim(tasktype),' task running on node=',trim(nodeName(rank))
100 format (a,i5,a,a,a,a)
  enddo
#endif
if (ios /= 0) then
  close (unitno)
endif
!SMS$SERIAL END

    return
  end subroutine print_nodeinfo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! printmem: Print max/min memory usage info across tasks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine printmem (str)
    character(len=*), intent(in) :: str                 ! string to print
    integer :: procsiz, rss, share, text, datastack     ! returned by gptlget_memusage
    integer :: rssmax, rssmin                           ! max/min across MPI ranks
    integer :: ret

#include <gptl.inc>

    ret = gptlget_memusage (procsiz, rss, share, text, datastack)
    rssmax = rss
    rssmin = rss
!SMS$REDUCE(rssmax,max)
!SMS$REDUCE(rssmin,min)
    write(6,*) str, ' Task memory usage info:'
    write(6,*)'Max,min process sizes=', rssmax, rssmin
    return
  end subroutine printmem
end module print_taskinfo
