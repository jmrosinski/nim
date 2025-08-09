integer function gptlpr_summary (comm)
  implicit none

  integer, intent(in) :: comm

  gptlpr_summary = 0
  return
end function gptlpr_summary

integer function gptlinitialize ()
  implicit none

  gptlinitialize = 0
  return
end function gptlinitialize

integer function gptlstart (name)
  implicit none

  character(len=*), intent(in) :: name

  gptlstart = 0
  return
end function gptlstart

integer function gptlstop (name)
  implicit none

  character(len=*), intent(in) :: name

  gptlstop = 0
  return
end function gptlstop

integer function gptlget_nregions (t, nregions)
  implicit none
  
  integer, intent(in) :: t
  integer, intent(out) :: nregions

  nregions = 0
  gptlget_nregions = 0
  return
end function gptlget_nregions

integer function gptlget_regionname (t, n, name)
  implicit none

  integer, intent(in) :: t
  integer, intent(in) :: n
  character(len=*), intent(out) :: name

  name = ' '
  gptlget_regionname = 0
  return
end function gptlget_regionname

integer function gptlget_wallclock (name, t, value)
  implicit none

  character(len=*), intent(in) :: name
  integer, intent(in) :: t
  real(8), intent(out) :: value

  value = 0.
  gptlget_wallclock = 0
  return
end function gptlget_wallclock

integer function gptlpr (mype)
  implicit none

  integer, intent(in) :: mype

  gptlpr = 0
  return
end function gptlpr

integer function gptlprint_memusage (string)
  implicit none

  character(len=*), intent(in) :: string
  gptlprint_memusage = 0
  return
end function gptlprint_memusage

integer function gptlget_memusage (size, rss, share, text, datastack)
  implicit none

  integer, intent(out) :: size, rss, share, text, datastack
  size = 0
  rss = 0
  share = 0
  text = 0
  datastack = 0
  gptlget_memusage = 0
  return
end function gptlget_memusage

integer function gptlprint_rusage (string)
  implicit none

  character(len=*), intent(in) :: string
  gptlprint_rusage = 0
  return
end function gptlprint_rusage
