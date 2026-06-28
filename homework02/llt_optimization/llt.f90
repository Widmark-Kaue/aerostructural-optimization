module llt

    ! Module that contains the necessary subroutines
    ! to apply the lifting-line theory (LLT) to a planar
    ! rectangular wing (constant chord) and variable
    ! twist along the span.
    !
    ! Aeronautics Institute of Technology
    ! Sao Jose dos Campos, SP, Brazil
    
    ! Author: Widmark Kauê Silva Cardoso
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
        integer :: ii
        real :: pi
        real :: panel_length, Sref
        real :: yc(nelem)                           ! Spanwise locations of the horseshoe vortices
        real :: wi(nelem), Veff(nelem)              ! Induced velocity at each horseshoe vortex and efective velocity at each wing section
        real :: alpha_ind(nelem), alpha_eff(nelem)  ! Induced angle of attack and effective angle of attack at each wing section
        real :: cli(nelem), lsi(nelem), lki(nelem)  ! Section lift coefficient, sectional lift force and Kutta-Joukowski force at each wing section
        real :: idx(nelem)                          ! Index array for the horseshoe vortices
        
        ! EXECUTION
        
        ! Constants
        pi = acos(-1.0)
        panel_length = span / real(nelem)
        Sref = span * chord
        idx = [(1.0 + (ii-1)*1.0, ii = 1, nelem)]

        ! ADD YOUR CODE HERE 
        yc  = 0.5*(panel_length - span) + (idx-1)*panel_length ! Spanwise locations of the horseshoe vortices

        ! Compute the induced velocity at each horseshoe vortex using the Biot-Savart law
        do ii = 1, nelem
            wi(ii) = panel_length/pi * sum(gama/(4*(yc(ii) - yc)**2 - panel_length**2)) ! Induced velocity at each horseshoe vortex
        end do

        ! Compute the induced AoA, effective AoA and effective velocity at each wing section
        alpha_ind = atan(wi/Vinf) ! Induced angle of attack at each wing section
        alpha_eff = alpha + twist + alpha_ind ! Effective angle of attack at each wing section
        Veff = sqrt(Vinf**2 + wi**2) ! Effective velocity at each wing section
        cli = cl0 + cla*alpha_eff ! Section lift coefficient at each wing section
        
        ! Compute lift force at each section 
        lsi = 0.5*rho_air*Veff**2*chord*cli ! Sectional lift force at each wing section
        lki = rho_air*Veff*gama ! Kutta-Joukowski force at each wing section

        ! Compute the LLT residuals at each wing section
        res_llt = lsi - lki ! LLT residuals at each wing section
        
        ! Compute the total lift and drag coefficients of the wing

        CL = sum(lki*panel_length*cos(alpha_ind)) / (0.5*rho_air*Vinf**2* Sref) ! Total lift coefficient of the wing
        CD = -sum(lki*panel_length*sin(alpha_ind)) / (0.5*rho_air*Vinf**2* Sref) ! Total drag coefficient of the wing

        ! END OF YOUR CODE

    end subroutine llt_main

end module llt
