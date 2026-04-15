! main.f90
! program that uses multiple derivative calculation methods

program main

  use mymodule, only : test_func, test_func_deriv, test_func_complex

  implicit none

  real :: x,h,f,fh
  real :: dfdx_an
  real :: dfdx_fd, error_fd
  complex :: fhc
  real :: dfdx_cs, error_cs
  ! Define test point
  x = -0.5

  ! Define step size
  h = 1e-8

  ! Compute analytical derivative
  call test_func_deriv(x,dfdx_an)
  print *,'dfdx  AN:',dfdx_an

  ! FINITE DIFFERENCE
  call test_func(x,f)
  call test_func(x+h, fh)

  dfdx_fd = (fh - f)/h
  error_fd = abs(1 - dfdx_fd/dfdx_an)
  print *,''
  print *,'dfdx  FD:',dfdx_fd
  print *,'error FD:',error_fd

  ! COMPLEX STEP
  call test_func_complex(complex(x, h), fhc)
  dfdx_cs = aimag(fhc)/h
  error_cs = abs(1 - dfdx_cs/dfdx_an)
  print *,''
  print *,'dfdx  CS:',dfdx_cs
  print *,'error CS:',error_cs

end program main
