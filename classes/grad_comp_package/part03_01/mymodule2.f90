! mymodule2.f90
! module with a multivariable function for derivative calculation

module mymodule2

  implicit none

contains

  subroutine test_func(x,y,z,f)

    implicit none
    real, intent(in) :: x,y,z
    real, intent(out) :: f
    real :: a

    a = y*cos(x) + z

    f = a**2 + 2.0/x

  end subroutine test_func

  !===================================

  subroutine test_func_complex(x,y,z,f)

    implicit none
    complex, intent(in) :: x,y,z
    complex, intent(out) :: f
    complex :: a

    a = y*cos(x) + z

    f = a**2 + 2.0/x

  end subroutine test_func_complex

  !===================================

  subroutine test_func_deriv(x,y,z,dfdx,dfdy,dfdz)

    ! Analytical derivative of the test function

    implicit none
    real, intent(in) :: x,y,z
    real, intent(out) :: dfdx,dfdy,dfdz

    dfdx = -y**2*sin(2.0*x) - 2.0*y*z*sin(x) - 2.0/x**2

    dfdy = 2*y*cos(x)**2 + 2.0*z*cos(x)

    dfdz = 2.0*y*cos(x) + 2.0*z

  end subroutine test_func_deriv

  !===================================

  subroutine test_func_md(x,xd,y,yd,z,zd,f,fd)

    implicit none
    real, intent(in) :: x,y,z
    real, intent(in) :: xd,yd,zd
    real, intent(out) :: f
    real, intent(out) :: fd
    real :: a
    real :: ad

    a = y*cos(x) + z

    ad = -y*sin(x)*xd + cos(x)*yd + zd

    f = a**2 + 2.0/x

    fd = 2.0*a*ad - 2.0/x**2*xd

  end subroutine test_func_md

  !===================================

  subroutine test_func_mb(x,xb,y,yb,z,zb,f,fb)

    implicit none
    real, intent(in) :: x,y,z
    real, intent(out) :: xb,yb,zb
    real, intent(out) :: f
    real, intent(in) :: fb
    real :: a
    real :: ab

    ! Forward pass
    a = y*cos(x) + z
    f = a**2 + 2.0/x

    ! Reverse pass

    ! Initialize seeds
    xb = 0.0
    yb = 0.0
    zb = 0.0
    ab = 0.0

    ! f = a**2 + 2.0/x
    ab = ab + 2*a*fb
    xb = xb - 2.0/x**2*fb

    ! a = y*cos(x) + z
    xb = xb - y*sin(x)*ab
    yb = yb + cos(x)*ab
    zb = zb + ab

  end subroutine test_func_mb

end module mymodule2
