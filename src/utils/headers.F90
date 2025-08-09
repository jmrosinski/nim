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

!These are the four routines that write, and verify the headers in the *.dat file.
!The I/O in all four routines must match.

subroutine WriteGlvlHeader(unit,glvl)
integer     ,intent(IN ) :: unit
integer     ,intent(IN ) :: glvl
character(16)            :: header
write(header,"('glvl =',I2)") glvl
write(unit) header
return
end subroutine WriteGlvlHeader

subroutine WriteCurveHeader(unit,curve)
integer     ,intent(IN ) :: unit
integer     ,intent(IN ) :: curve
character(16)            :: header
write(header,"('curve =',I2)") curve
write(unit) header
return
end subroutine WriteCurveHeader

subroutine TestGlvlHeader(unit,FileName,RoutineName,glvl)
integer     ,intent(IN) :: unit
character(*),intent(IN) :: FileName
character(*),intent(IN) :: RoutineName
integer     ,intent(IN) :: glvl
integer                 :: glvlHeader
character(16)           :: header
character(80)           :: FMT
FMT ="Error in ',a,' unit ',i0,' glvl=,i0,' does not match header glvl=',i0"
read(unit) header
read(header,"(6x,I2)") glvlHeader
if(glvl /= glvlHeader) then
  write(6,FMT) RoutineName,unit,glvl,glvlHeader
endif
end subroutine TestGlvlHeader

subroutine TestCurveHeader(unit,FileName,RoutineName,curve)
integer     ,intent(IN) :: unit
character(*),intent(IN) :: FileName
character(*),intent(IN) :: RoutineName
integer     ,intent(IN) :: curve
integer                 :: curveHeader
character(16)           :: header
character(80)           :: FMT
FMT ="Error in ',a,' unit ',i0,' curve=,i0,' does not match header curve=',i0"
read(unit) header
read(header,"(7x,I2)") curveHeader
if(curve /= curveHeader) then
  write(6,FMT) RoutineName,unit,curve,curveHeader
endif
end subroutine TestCurveHeader

