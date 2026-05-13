! mymodule.f90
! module with function definitions for derivative calculation

module mymodule

  implicit none

contains

  subroutine test_func(x,f)

    implicit none
    real, intent(in) :: x
    real, intent(out) :: f

    if (x > 0.0) then
       f = x**2/(1+x)
    else
       f = abs(x)**3
    end if

  end subroutine test_func

  !===================================

   subroutine test_func_complex(x,f)

    implicit none
    complex, intent(in) :: x
    complex, intent(out) :: f

    if (real(x) > 0.0) then
       f = x**2/(1+x)
    else
       f = abs_cs(x)**3
    end if

  end subroutine test_func_complex

  !===================================

  subroutine test_func_deriv(x,dfdx)

    ! Analytical derivative of the test function

    implicit none
    real, intent(in) :: x
    real, intent(out) :: dfdx

    if (x > 0.0) then
       dfdx = 2*x/(1+x) - x**2/(1+x)**2
    else
       dfdx = -3*x**2
    end if

  end subroutine test_func_deriv

  
  !===================================

  function abs_cs(z)
    ! No function não precisa dizer intent(in) ou intent(out)
    ! A function identifica a variável de saída pelo nome da função, que é abs_cs nesse caso

    implicit none
    complex :: z
    complex :: abs_cs


    if (real(z) > 0.0) then
       abs_cs = z
    else
       abs_cs = -z
    end if

  end function abs_cs


end module mymodule
