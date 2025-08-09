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

module namelistdata
!******************************************************************************
! NIM namelist data that must be read by Fortran programs and also used by 
! scripts is read by routines defined here.  
!******************************************************************************
  use ReadNamelist
  implicit none

contains

  subroutine GetNprocs  (nprocs, ret)
    integer, intent(out) :: nprocs
    integer, intent(out) :: ret
    integer :: rc
    integer :: ComputeTasks
    character(len=13) :: nlfile='./NIMnamelist'
    namelist /COMPUTETASKnamelist/ ComputeTasks
    open (unit=10,file=nlfile,status='old',action='read',iostat=rc)
    if (rc /= 0) then
      write (6,'(a,a)') 'GetNprocs: ERROR opening ',trim(nlfile)
      ret = rc
      return
    endif
    read (unit=10, NML=COMPUTETASKnamelist, iostat=rc)
    if (rc /= 0) then
      write(6,'(a,a)') 'GetNprocs: ERROR reading COMPUTETASKnamelist from ', &
                          trim(nlfile)
      ret = rc
      return
    endif
    close (10)
    nprocs = ComputeTasks
    ret = 0
    return
  end subroutine GetNprocs

  subroutine GetMaxQueueTime (QueueTime, ret)
    character(len=8), intent(out) :: QueueTime
    integer, intent(out) :: ret
    integer :: rc
    call readnl (rc)
    if (rc /= 0) then
      QueueTime = 'XXXXXXXX'
      ret = -1
      print*, 'GetMaxQueueTime: bad return from readnl'
      return
    end if
    QueueTime = MaxQueueTime
    ret = 0
    return
  end subroutine GetMaxQueueTime

  subroutine ReturnGLVL (glvlout, ret)
    integer, intent(out) :: glvlout
    integer, intent(out) :: ret
    integer :: rc
    call readnl (rc)
    if (rc /= 0) then
      glvlout = -1
      ret = -1
      print*, 'ReturnGLVL: bad return from readnl'
      return
    end if
    glvlout = glvl
    ret = 0
    return
  end subroutine ReturnGLVL

  subroutine ReturnNIP (nipout, ret)
    integer,intent(OUT) :: nipout
    integer ret, i

    integer :: rc

    call readnl (rc)
    if (rc /= 0) then
      nipout =  -1
      ret = -1
      write(6,*)'ReturnNIP: bad return from readnl'
      return
    end if

    nipout = 1
    do i = 1, glvl
        nipout = nipout * SubdivNum(i)
    enddo
    nipout = 10 * nipout * nipout + 2
  end subroutine ReturnNIP

  subroutine ReturnNZ (nzout, ret)
    integer, intent(out) :: nzout
    integer, intent(out) :: ret
    integer :: rc
    call readnl (rc)
    if (rc /= 0) then
      nzout = -1
      ret = -1
      print*, 'ReturnNZ: bad return from readnl'
      return
    end if
    nzout = nz
    ret = 0
    return
  end subroutine ReturnNZ

  subroutine ReturnGridType (gtout, ret)
    integer,intent(OUT) :: gtout
    integer,intent(OUT) :: ret
    integer :: rc
    call readnl(rc)
    if (rc /= 0) then
      gtout = -1
      ret = -1
      print*, 'ReturnGridType: bad return from readnl'
      return
    end if
    gtout = gtype
    ret = 0
    return
  end subroutine ReturnGridType

  subroutine ReturnYYYYmmddhhmm(yyyymmddhhmmout, ret)
    character(12), intent(out) :: yyyymmddhhmmout
    integer, intent(out) :: ret

    integer :: rc

    call readnl (rc)
    if (rc /= 0) then
      yyyymmddhhmmout = ''
      ret = -1
      write(6,*)'ReturnYYYYmmddhhmm: bad return from readnl'
      return
    end if

    yyyymmddhhmmout = yyyymmddhhmm
  end subroutine ReturnYYYYmmddhhmm
  
  subroutine GetDataDir(datadirout, ret)
    character(128), intent(out) :: datadirout
    integer, intent(out) :: ret

    integer :: rc

    ret = 0

    call readnl (rc)
    if (rc /= 0) then
      datadirout = ''
      ret = -1
      write(6,*)'GetDataDir: bad return from readnl'
      return
    end if

    datadirout = datadir
  end subroutine GetDataDir

  subroutine GetDataGrd(datagrdout, ret)
    character(128), intent(out) :: datagrdout
    integer, intent(out) :: ret

    integer :: rc

    ret = 0

    call readnl (rc)
    if (rc /= 0) then
      datagrdout = ''
      ret = -1
      write(6,*)'GetDataGrd: bad return from readnl'
      return
    end if

    datagrdout = datadir
  end subroutine GetDataGrd

  subroutine ReturnPhysics (physicsOut, ret)
    character(12),intent(out) :: physicsOut ! Physical package: GRIMS or GFS or none for no physics
    integer      ,intent(out) :: ret
    integer :: rc
    call readnl (rc)
    if (rc /= 0) then
      physicsOut = 'Error'
      ret = -1
      print*, 'ReturnPhysics: bad return from readnl'
      return
    end if
    physicsOut = physics
    ret = 0
    return
  end subroutine ReturnPhysics

  subroutine GetRestartDir(RestartDirOut, ret)
    character(128), intent(out) :: RestartDirOut
    integer, intent(out) :: ret

    integer :: rc

    ret = 0

    call readnl (rc)
    if (rc /= 0) then
      RestartDirOut = ''
      ret = -1
      write(6,*)'GetRestartDir: bad return from readnl'
      return
    end if
    RestartDirOut = RestartDir

  end subroutine GetRestartDir

  subroutine ReturnRESTARTBEGIN (restartbeginOut, ret)
    integer, intent(out) :: restartbeginOut
    integer, intent(out) :: ret
    integer :: rc
    call readnl (rc)
    if (rc /= 0) then
      restartbeginOut = -1
      ret = -1
      print*, 'ReturnRESTARTBEGIN: bad return from readnl'
      return
    end if
    restartbeginOut = RestartBegin
    ret = 0
    return
  end subroutine ReturnRESTARTBEGIN

end module namelistdata
