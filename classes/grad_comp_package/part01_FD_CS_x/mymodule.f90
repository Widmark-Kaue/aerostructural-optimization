! mymodule.f90
! module with function definitions for derivative calculation

module mymodule

  implicit none

contains

  subroutine test_func(x,f)

    implicit none
    real, intent(in) :: x
    real, intent(out) :: f

    f = x**2/(1+x)

  end subroutine test_func

  !===================================

  subroutine test_func_complex(x,f)

    implicit none
    complex, intent(in) :: x
    complex, intent(out) :: f

    f = x**2/(1+x)

  end subroutine test_func_complex

  !===================================

  subroutine test_func_deriv(x,dfdx)

    ! Analytical derivative of the test function

    implicit none
    real, intent(in) :: x
    real, intent(out) :: dfdx

    dfdx = 2*x/(1+x) - x**2/(1+x)**2

  end subroutine test_func_deriv

end module mymodule
