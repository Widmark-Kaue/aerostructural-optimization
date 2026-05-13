! fem_module.f90
! module with FEM problem implementation

module fem_module

  implicit none

contains

  subroutine get_residuals(n_bars, A, d, res)
    ! This subroutine computes the residuals of the FEM problem
    ! of the bar with uniform load
    !
    ! INPUTS
    ! n_bars: integer -> Number of bars in the discretization
    ! A : real(n_bars) -> Areas of each bar
    ! d : real(n_bars) -> Tip displacement of each bar
    !
    ! OUTPUTS
    ! res : real(n_bars) -> Residuals of the FEM problem

    implicit none

    ! Input variables
    integer, intent(in) :: n_bars
    real, intent(in) :: A(n_bars), d(n_bars)

    ! Output variables
    real, intent(out) :: res(n_bars)

    ! Working variables
    integer :: ii
    real :: K(n_bars, n_bars), f(n_bars), Kd(n_bars)

    ! EXECUTION

    !!! Build stiffness matrix
    K = 0.0

    ! Contribution of the first bar
    K(1,1) = A(1)

    ! Contribution of remaining bars
    do ii = 2,n_bars
       K(ii-1,ii-1) = K(ii-1,ii-1) + A(ii)
       K(ii-1,ii) = K(ii-1,ii) - A(ii)
       K(ii,ii-1) = K(ii,ii-1) - A(ii)
       K(ii,ii) = K(ii,ii) + A(ii)
    end do

    ! Multiply by the 1/L_bar factor
    ! In this problem L_bar = 1/n_bars
    K = K*n_bars

    !!! Build force vector (load/length = 1)

    f = 1.0/n_bars
    f(n_bars) = 0.5/n_bars

    !!! Compute residuals

    call matdotvec(n_bars, K, d, Kd)

    res = Kd - f

  end subroutine get_residuals

  !===================================

  subroutine get_functions(n_bars, A, d, C, V)
    ! This subroutine computes the functions of interest of the FEM problem
    ! of the bar with uniform load
    !
    ! INPUTS
    ! n_bars: integer -> Number of bars in the discretization
    ! A : real(n_bars) -> Areas of each bar
    ! d : real(n_bars) -> Tip displacement of each bar
    !
    ! OUTPUTS
    ! C : real -> Compliance (sum of f*d)
    ! V : real -> Bar volume

    implicit none

    ! Input variables
    integer, intent(in) :: n_bars
    real, intent(in) :: A(n_bars), d(n_bars)

    ! Output variables
    real, intent(out) :: C, V

    ! Working variables
    real :: f(n_bars)
    integer :: ii

    ! EXECUTION

    !!! Build force vector (load/length = 1)

    f = 1.0/n_bars
    f(n_bars) = 0.5/n_bars

    !!! Compliance
    C = 0

    !!! Volume
    V = 0

    do ii = 1,n_bars

      C = C + f(ii)*d(ii)

      V = V + A(ii)

    end do

    V = V/n_bars

  end subroutine get_functions

  !===================================

  subroutine matdotvec(n, K, d, Kd)
    ! Thsi subroutines does the dot product between matrix K
    ! and vector d
    !
    ! INPUTS
    ! n : integer -> matrix size
    ! K : real(n,n) -> matrix
    ! d : real(n) -> vector
    !
    ! OUTPUTS
    ! Kd : real(n) -> K*d


    implicit none

    ! Input variables
    integer, intent(in) :: n
    real, intent(in) :: K(n,n), d(n)

    ! Output variables
    real, intent(out) :: Kd(n)

    ! Working variables
    integer :: ii

    ! EXECUTION
    Kd = 0.0
    do ii = 1,n
       Kd = Kd + K(:,ii)*d(ii)
    end do

  end subroutine matdotvec

  !===================================

  subroutine main(n_bars, A, d, res, C, V)
    ! This is the main subroutine that calls two important subroutines
    ! to force their differentiation in a single Tapenade call,
    ! since the Tapenade web interface accepts a single
    ! subroutine name at a time.
    !
    ! INPUTS
    ! n_bars: integer -> Number of bars in the discretization
    ! A : real(n_bars) -> Areas of each bar
    ! d : real(n_bars) -> Tip displacement of each bar
    !
    ! OUTPUTS
    ! res : real(n_bars) -> Residuals of the FEM problem
    ! C : real -> Compliance (sum of f*d)
    ! V : real -> Bar volume

    implicit none

    ! Input variables
    integer, intent(in) :: n_bars
    real, intent(in) :: A(n_bars), d(n_bars)

    ! Output variables
    real, intent(out) :: res(n_bars)
    real, intent(out) :: C, V

    ! EXECUTION

    call get_residuals(n_bars, A, d, res)

    call get_functions(n_bars, A, d, C, V)

  end subroutine main

end module fem_module
