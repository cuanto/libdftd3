! Simple demonstration of how to use dftd3 as a library.
!
! Copyright (C) 2016, Bálint Aradi
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 1, or (at your option)
! any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! For the GNU General Public License, see <http://www.gnu.org/licenses/>
!

! Tests the dftd3 API by calculating the dispersion for a small DNA fragment.
!
program testapi2

  use dftd3_api
  implicit none

  ! Working precision
  integer, parameter :: wp = kind(1.0d0)

  ! Same conversion factor as in dftd3
  real(wp), parameter :: AA__Bohr = 1.0_wp / 0.52917726_wp

  integer, parameter :: nAtoms = 12
  
  ! Coordinates in Angstrom 
  ! They must be converted to Bohr before passed to dftd3
  real(wp), parameter :: coords(3, nAtoms) = reshape([ &
 &  -0.471925, -0.471925, -1.859111, &
 &   0.471925,  0.471925, -1.859111, &
 &  -0.872422, -0.872422, -0.936125, &
 &   0.872422,  0.872422, -0.936125, &
 &  -0.870464, -0.870464, -2.783308, &
 &   0.870464,  0.870464, -2.783308, &
 &  -0.471925,  0.471925,  1.859111, &
 &   0.471925, -0.471925,  1.859111, &
 &  -0.872422,  0.872422,  0.936125, &
 &   0.872422, -0.872422,  0.936125, &
 &  -0.870464,  0.870464,  2.783308, &
 &   0.870464, -0.870464,  2.783308  &
       &] * AA__Bohr, [3, nAtoms])
  
  integer, parameter :: species(nAtoms) = [1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2]

  integer, parameter :: nSpecies = 2
  character(2), parameter :: speciesNames(nSpecies) = [ 'C ', 'H ']
  
  type(dftd3_input) :: input
  type(dftd3_calc) :: dftd3
  integer :: atnum(nAtoms), i
  real(wp) :: edisp
  real(wp) :: grads(3, nAtoms)

  ! Initialize dftd3
  call dftd3_init(dftd3, input)

  ! Choose functional. Alternatively you could set the parameters manually
  ! by the dftd3_set_params() function.
  call dftd3_set_functional(dftd3, func='hf', version=3, tz=.false.)

  ! Convert species name to atomic number for each atom
  atnum(:) = get_atomic_number(speciesNames(species))

  ! Calculate dispersion and gradients for non-periodic case
  call dftd3_dispersion(dftd3, coords, atnum, edisp, grads)
  write(*, "(A)") "*** Dispersion for non-periodic case"
  write(*, "(A,ES20.12)") "Energy [au]:", edisp
  write(*, "(A)") "Gradients [au]:"
  write(*, "(3ES20.12)") grads
  write(*, *)

end program testapi2
