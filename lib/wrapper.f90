subroutine wrapper(natoms, coords, itype, fun, version, tz, edisp, grads) bind(c)

  use iso_c_binding, only: c_int, c_double, c_null_char, c_char
  use dftd3_api, only: dftd3_calc, dftd3_input, dftd3_init, dftd3_dispersion, &
                       dftd3_set_functional
  implicit none

  integer(kind=c_int), intent(in), value :: natoms
  real(kind=c_double) :: coords(natoms,3)
  integer(kind=c_int), intent(in) :: itype(natoms)
  character(kind=c_char,len=1), intent(in) :: fun(*)
  integer(kind=c_int), intent(in), value :: version
  integer(kind=c_int), intent(in), value :: tz
  real(kind=c_double), intent(out) :: edisp
  real(kind=c_double), intent(out) :: grads(3,natoms)

  type(dftd3_calc) :: dftd3
  type(dftd3_input) :: input
  real(kind=c_double) :: coordsr(3,natoms)
  logical :: ltz
  character(len=:), allocatable :: func
  integer(kind=c_int) :: i, nchars

  i = 1
  do
    if (fun(i) == c_null_char) exit
    i = i + 1
  end do
  nchars = i - 1
  allocate(character(len=nchars) :: func)
  func = transfer(fun(1:nchars), func)

  do i = 1,natoms
    coordsr(:,i) = coords(i,:)
  end do 

  if (tz.eq.0) ltz = .false.
  if (tz.eq.1) ltz = .true.

  ! Initialize dftd3
  call dftd3_init(dftd3, input)

  ! Choose functional.
  call dftd3_set_functional(dftd3, func=func, version=version, tz=ltz)

  ! Calculate dispersion and gradients for non-periodic case
  call dftd3_dispersion(dftd3, coordsr, itype, edisp, grads)

end subroutine wrapper
