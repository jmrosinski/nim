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

module module_header

  use module_control,only:nip,glvl,yyyymmddhhmm
  use ReadNamelist, only: nz
!  use namelistdata,only:glvl,nz,yyyymmddhhmm

  implicit none
  public header
  private

! Parameters

  integer,parameter::header_cols=80
  integer,parameter::header_rows=10
  integer,parameter::header_len=header_cols*header_rows

! Module variables

  character(len=header_cols)::line
  character(len=header_len)::h
  integer::pos

contains

  function header(varname,levels,its,time,timeunit)
    implicit none
    character*(*),intent(in)::varname,timeunit
    character(len=header_len)::header
    integer,intent(in)::its,levels,time
    pos=1
    h=' '
    write (line,'(a,a,a,a)') 'NIM ',varname,' Forecast initial time YYYYMMDDHHMM: ',yyyymmddhhmm
    call append
!KY    write (line,'(a,i0,a,i0,a,i0,a,i0,a,a)') 'Level ',nz,', GLVL= ',glvl,', Step ',its,', ',time,' ','hours'
    write (line,'(a,i0,a,i0,a,a,a,i0,a,i0,a,a)') 'Level ',nz,', GLVL= ',glvl,', Time Unit = ',timeunit,', Step ',its,', ',time,' ','hours'
    call append
    write (line,'(a,i0,a,i0)') 'dim1= ',levels,', nip= ',nip
    call append
    write (line,'(i0)') 4
    call append
    write (line,'(i0)') 5
    call append
    write (line,'(i0)') 6
    call append
    write (line,'(i0)') 7
    call append
    write (line,'(i0)') 8
    call append
    write (line,'(i0)') 9
    call append
    write (line,'(i0)') 10
    call append
    header=h
  end function header

  subroutine append
    implicit none
    integer::i,j
    if (pos.ge.header_len) then
      write (*,'(a)') 'ERROR: Attempt to write past end of header.'
      stop
    endif
    do i=1,len(line)
      j=pos+i-1
      h(j:j)=line(i:i)
    enddo
    pos=(int((pos+header_cols)/header_cols)*header_cols)+1
  end subroutine append

end module module_header
