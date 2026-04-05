module llt

    ! Module that contains the necessary subroutines
    ! to apply the lifting-line theory (LLT) to a planar
    ! rectangular wing (constant chord) and variable
    ! twist along the span.
    !
    ! Aeronautics Institute of Technology
    ! Sao Jose dos Campos, SP, Brazil
    
    implicit none
    
contains

    subroutine llt_main(nelem, twist, gama, span, &
    chord, cl0, cla, alpha, Vinf, rho_air, res_llt, CL, CD)
    
    ! This subroutine receives the rectangular wing geometry and flight conditions
    ! to compute residuals and aerodynamic forces according to the LLT method.
    !
    ! INPUTS
    !
    ! nelem: integer -> Number of horseshoe vortices
    ! twist: real(nelem) -> Incidence angle of each wing section [rad]
    ! gama: real(nelem) -> Circulation of each horseshoe vortex [m2/s]
    ! span: real -> Total wing span [m]
    ! chord: real -> Wing chord [m]
    ! cl0: real -> Airfoil lift coefficient at zero angle of attack
    ! cla: real -> Airfoil lift coefficient slope [1/rad]
    ! alpha: real -> Wing angle of attack [rad]
    ! Vinf: real -> Freestream airspeed [m/s]
    ! rho_air: real -> Freestream density [kg/m3]
    !
    ! OUTPUTS
    !
    ! res_llt: real(nelem) -> LLT residuals at each wing section [N/m]
    ! CL: real -> Total lift coefficient of the wing
    ! CD: real -> Total drag coefficient of the wing
    
    implicit none
    
    ! Input variables
    integer, intent(in) :: nelem
    real, intent(in) :: twist(nelem), gama(nelem)
    real, intent(in) :: span, chord, cl0, cla, alpha, Vinf, rho_air
    
    ! Output variables
    real, intent(out) :: res_llt(nelem)
    real, intent(out) :: CL, CD
    
    ! Working variables
    real :: pi
    ! EXECUTION
    
    ! Constants
    pi = acos(-1.0)
    
    ! ADD YOUR CODE HERE 

    end subroutine llt_main
    
end module llt
