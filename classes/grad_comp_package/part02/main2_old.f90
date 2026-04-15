! main2.f90
! program that uses multiple derivative calculation methods

program main2

  use mymodule2, only : test_func, test_func_deriv, test_func_complex

  implicit none

  real :: x,y,z
  real :: h,f,fh
  real :: dfdx_an, dfdy_an, dfdz_an
  real :: dfdx_fd, errorx_fd, dfdy_fd, errory_fd, dfdz_fd, errorz_fd
  complex :: fhc
  real :: dfdx_cs, errorx_cs, dfdy_cs, errory_cs, dfdz_cs, errorz_cs

  ! Define test point
  x = 1.5
  y = 2.5
  z = 10.0

  ! Define step size
  h = 1e-6

  ! Compute analytical derivative
  call test_func_deriv(x,y,z,dfdx_an,dfdy_an,dfdz_an)
  print *,'dfdx  AN:',dfdx_an,dfdy_an,dfdz_an

  ! FINITE DIFFERENCE
  call test_func(x, y, z, f)

  ! dfdx
  ???
  dfdx_fd = ???
  errorx_fd = abs(1 - dfdx_fd/dfdx_an)

  ! dfdy
  ???
  dfdy_fd = ???
  errory_fd = abs(1 - dfdy_fd/dfdy_an)

  ! dfdz
  ???
  dfdz_fd = ???
  errorz_fd = abs(1 - dfdz_fd/dfdz_an)

  print *,''
  print *,'dfdx  FD:',dfdx_fd,dfdy_fd,dfdz_fd
  print *,'error FD:',errorx_fd,errory_fd,errorz_fd

  ! COMPLEX-STEP

  ! dfdx
  call test_func_complex(complex(x,h), complex(y,0.0), complex(z,0.0), fhc)
  dfdx_cs = imag(fhc)/h
  errorx_cs = abs(1 - dfdx_cs/dfdx_an)

  ! dfdy
  ???
  dfdy_cs = imag(fhc)/h
  errory_cs = abs(1 - dfdy_cs/dfdy_an)

  ! dfdz
  ???
  dfdz_cs = imag(fhc)/h
  errorz_cs = abs(1 - dfdz_cs/dfdz_an)

  print *,''
  print *,'dfdx  CS:',dfdx_cs,dfdy_cs,dfdz_cs
  print *,'error CS:',errorx_cs,errory_cs,errorz_cs

  ! FORWARD AD

  ! dfdx
  !xd = ???
  !yd = ???
  !zd = ???
  !call test_func_md(x,xd,y,yd,z,zd,f,fd)
  !dfdx_md = ???
  !errorx_md = abs(1 - dfdx_md/dfdx_an)

  ! dfdy
  !xd = ???
  !yd = ???
  !zd = ???
  !call test_func_md(x,xd,y,yd,z,zd,f,fd)
  !dfdy_md = ???
  !errory_md = abs(1 - dfdy_md/dfdy_an)

  ! dfdz
  !xd = ???
  !yd = ???
  !zd = ???
  !call test_func_md(x,xd,y,yd,z,zd,f,fd)
  !dfdz_md = ???
  !errorz_md = abs(1 - dfdz_md/dfdz_an)

  !print *,''
  !print *,'dfdx  MD:',dfdx_md,dfdy_md,dfdz_md
  !print *,'error MD:',errorx_md,errory_md,errorz_md

  ! REVERSE AD

  ! Initialize reverse seeds
  !fb = ???
  !call test_func_mb(x,xb,y,yb,z,zb,f,fb)
  !dfdx_mb = ???
  !dfdy_mb = ???
  !dfdz_mb = ???
  !errorx_mb = abs(1 - dfdx_mb/dfdx_an)
  !errory_mb = abs(1 - dfdy_mb/dfdy_an)
  !errorz_mb = abs(1 - dfdz_mb/dfdz_an)

  !print *,''
  !print *,'dfdx  MB:',dfdx_mb,dfdy_mb,dfdz_mb
  !print *,'error MB:',errorx_mb,errory_mb,errorz_mb

end program main2
