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

end module mymodule2
