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

!*********************************************************************
! Apply a perturbation to an input field, using a rectangular distribution 
! bounded by "pertlim". Results will be invariant across processor counts
!
!   J. Rosinski (Nov, 2011)
!*********************************************************************

subroutine perturb (var, pertlim)
  use module_control, only: nip, ims, ime, ips, ipe
  use ReadNamelist, only: nz
  use kinds, only: rt
  implicit none

  real(rt), intent(in) :: pertlim        ! perturbation magnitude

  real(rt), intent(inout) :: var(nz,ims:ime) ! array to be perturbed
  real(rt) :: val(nip)              ! random numbers
  integer :: ipn
  integer :: ret

  call random_number (val(:))     ! generate random numbers
!$OMP PARALLEL DO
  do ipn=ips,ipe
    var(:,ipn) = var(:,ipn)*(1._rt + pertlim*(0.5_rt - val(ipn)))  ! apply the perturbation
  end do
!$OMP END PARALLEL DO

  return
end subroutine perturb
