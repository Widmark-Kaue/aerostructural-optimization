module llt_module

  ! Module that holds necessary subroutines for
  ! the generalized lifting line method.
  ! The formulation is based on:
  !   Modern Adaptation of Prandtl's Classic Lifting-Line Theory
  !   W. F. Phillips and D. O. Snyder
  !   Journal of Aircraft Vol. 37, No. 4, July-August 2000

  implicit none

  real :: zero = 0.0
  real :: pi = 3.1415926535897932384626
  real :: eps = epsilon(zero)

contains

  subroutine induced_vel(Uinf, X1, X2, Xp, gamma_ind)

    ! This subroutine computes the velocity induced by
    ! a horseshoe vortex with kinks at X1 and X2 at the
    ! point X.
    !
    ! INPUTS
    !
    ! Uinf: real(3) -> Unitary vector along the free-stream direction
    ! X1: real(3) -> Coordinates of the first kink of the horseshoe vortex
    ! X2: real(3) -> Coordinates of the second kink of the horseshoe vortex
    ! X: real(3) -> Coordinates where we want the induced velocity
    !
    ! OUTPUTS
    !
    ! gamma_ind: real(3) -> Aerodynamic influence vector induced
    !                       by the horseshoe vortex

    implicit none

    ! Input variables
    real, intent(in) :: Uinf(3), X1(3), X2(3), Xp(3)

    ! Output variables
    real, intent(out) :: gamma_ind(3)

    ! Working variables
    real :: r1(3), r2(3), V1(3), V2(3), V3(3)
    real :: UinfCr1(3), UinfCr2(3), r1Cr2(3)
    real :: r1m, r2m, UinfDr1, UinfDr2, r1Dr2, den

    ! EXECUTION

    ! Compute vectors between the kinks and the point of interest
    r1 = Xp - X1
    r2 = Xp - X2

    ! Compute the magnitude of these vectors
    call norm(r1, r1m)
    call norm(r2, r2m)

    ! Compute cross products
    call cross(Uinf, r1, UinfCr1)
    call cross(Uinf, r2, UinfCr2)
    call cross(r1, r2, r1Cr2)

    ! Compute dot products
    call dot(Uinf, r1, UinfDr1)
    call dot(Uinf, r2, UinfDr2)
    call dot(r1, r2, r1Dr2)

    ! Compute contributions of each vortex segment
    V1 = UinfCr2/r2m/(r2m - UinfDr2)

    den = r1m*r2m + r1Dr2
    if (abs(den) > eps) then
       V2 = (r1m+r2m)*r1Cr2/r1m/r2m/den
    else
       ! Bound vortex induces no velocity along its length
       V2 = zero
    end if

    V3 = -UinfCr1/r1m/(r1m - UinfDr1)

    ! Compute overall induced velocity (without the circulation: V/gamma)
    gamma_ind = 0.25*(V1 + V2 + V3)/pi

  end subroutine induced_vel

  !============================================================

  subroutine get_AIC_matrix(n_vort, X, Uinf, AIC)

    ! This subroutine computes the Aerodynamic Influence Coefficients (AIC)
    ! matrix that gives the adimensionalized induced velocities in
    ! all bound segments in X, that is: Vind = AIC*Gamma.
    ! Vind(1:3,ii) = AIC(1:3,:,ii)*Gamma(:)
    !
    ! INPUTS
    !
    ! n_vort: integer -> Number of horseshoe vortices.
    ! X: real(3,n_vort+1) -> Coordinates of the kinks of the horseshoe vortices.
    ! Uinf: real(3) -> Unitary vector along the free-stream direction.
    !
    ! OUTPUTS
    !
    ! AIC: real(3,n_vort,n_vort) -> Aerodynamic Influence matrix.

    implicit none

    ! Input variables
    integer, intent(in) :: n_vort
    real, intent(in) :: X(3,n_vort+1)
    real, intent(in) :: Uinf(3)

    ! Output variables
    real, intent(out) :: AIC(3,n_vort,n_vort)

    ! Working variables
    integer :: ii, jj
    real :: X1(3), X2(3), Xcen(3), gamma_ind(3)

    ! EXECUTION

    ! Loop over all vortices
    do ii=1,n_vort

       ! Get kink coordinates of the current vortex
       X1 = X(:,ii)
       X2 = X(:,ii+1)

       ! Compute the center of the bound vortex
       Xcen = 0.5*(X1 + X2)

       ! Loop over all bound vortices to compute the influence coefficients
       do jj=1,n_vort

          ! Get kink coordinates of the current vortex
          X1 = X(:,jj)
          X2 = X(:,jj+1)

          ! The induced velocity, normalized by gamma, is the influence
          ! coefficient.
          call induced_vel(Uinf, X1, X2, Xcen, gamma_ind)
          AIC(:,jj,ii) = gamma_ind

       end do

    end do

  end subroutine get_AIC_matrix

  !============================================================

  subroutine get_geom_vectors(n_vort, X, twist, chords, Uai, Uni, si, areas)

    ! This subroutine computes some geometric properties of the
    ! aerodynamic strips.
    !
    ! INPUTS
    !
    ! n_vort: integer -> Number of horseshoe vortices.
    ! X: real(3,n_vort+1) -> Coordinates of the kinks of the horseshoe vortices.
    ! twist: real(n_vort) -> Array of local incidence angles.
    ! chords: real(n_vort) -> Array of local chords.
    !
    ! OUTPUTS
    !
    ! Uai: real(3,n_vort) -> Unitary vectors along local chords.
    ! Uni: real(3,n_vort) -> Unitary vectors along local normals.
    ! si: real(3,n_vort) -> Vectors along bound vortices. They are not unitary.
    ! areas: real(n_vort) -> Areas of each panel.

    implicit none

    ! Input variables
    integer, intent(in) :: n_vort
    real, intent(in) :: X(3,n_vort+1)
    real, intent(in) :: twist(n_vort), chords(n_vort)

    ! Output variables
    real, intent(out) :: Uai(3,n_vort), Uni(3,n_vort)
    real, intent(out) :: si(3,n_vort), areas(n_vort)

    ! Working variables
    integer :: ii
    real :: X1(3), X2(3), Usi(3), ni(3)
    real :: nim, sim

    ! EXECUTION

    ! Loop over all vortices
    do ii=1,n_vort

       ! Get kink coordinates of the current vortex
       X1 = X(:,ii)
       X2 = X(:,ii+1)

       ! Compute the vector along the local span
       si(:,ii) = X2 - X1
       call normalize(si(:,ii), Usi)

       ! Create a unitary vector along the chordline
       ! pointing to the trailing edge
       Uai(1,ii) = cos(twist(ii))
       Uai(2,ii) = zero
       Uai(3,ii) = -sin(twist(ii))

       ! Get the normal vector as the normal between the span
       ! and the chord vector
       call cross(Uai(:,ii), Usi, ni)
       call normalize(ni, Uni(:,ii))

       ! Compute the local area of the panel.
       ! The norm of the normal vector is the sine of the angle.
       call norm(ni, nim)
       call norm(si(:,ii), sim)
       areas(ii) = nim*chords(ii)*sim

    end do

  end subroutine get_geom_vectors

  !============================================================

  subroutine get_local_vels(n_vort, Gama, AIC, Vinf, Uai, Uni, Vlocal, alphaLocal)

    ! This subroutine computes the local velocity an angle of attack at each
    ! bound vortex for the given vorticity values.

    implicit none

    ! Input variables
    integer, intent(in)  :: n_vort
    real, intent(in) :: Gama(n_vort), AIC(3,n_vort,n_vort)
    real, intent(in) :: Vinf(3), Uai(3,n_vort), Uni(3,n_vort)

    ! Output variables
    real, intent(out) :: Vlocal(3,n_vort), alphaLocal(n_vort)

    ! Working variables
    integer :: ii,jj
    real :: num, den

    ! EXECUTION

    ! Do the matrix operation
    do ii=1,n_vort

       ! Start with the free-strem contribution
       Vlocal(:,ii) = Vinf

       ! Now add the velocities induced by each vortex
       do jj=1,n_vort
          Vlocal(:,ii) = Vlocal(:,ii) + AIC(:,jj,ii)*Gama(jj)
       end do

       ! Compute the local angle of attack
       call dot(Vlocal(:,ii), Uni(:,ii), num)
       call dot(Vlocal(:,ii), Uai(:,ii), den)
       alphaLocal(ii) = atan(num/den)

    end do

  end subroutine get_local_vels

  !============================================================

  subroutine compute_circulation_forces(n_vort, rho, Vlocal, si, Gama, Fcirc)

    ! This subroutine uses the circulation values to compute aerodynamic
    ! forces on every bound vortex.
    !
    ! INPUTS

    implicit none

    ! Input variables
    integer, intent(in) :: n_vort
    real, intent(in) :: rho, Vlocal(3,n_vort), si(3,n_vort)
    real, intent(in) :: Gama(n_vort)

    ! Output variables
    real, intent(out) :: Fcirc(3,n_vort)

    ! Working variables
    integer :: ii
    real :: li(3)

    ! EXECUTION

    do ii = 1,n_vort

       ! L=rho*V*gamma*l
       call cross(Vlocal(:,ii), si(:,ii), li)
       Fcirc(:,ii) = rho*Gama(ii)*li

    end do

  end subroutine compute_circulation_forces

  !============================================================

  subroutine compute_airfoil_forces(n_vort, rho, cl0, cla, areas, Vlocal, alphaLocal, si, Fairf)

    ! This subroutine compares the 2D forces and the forces given by
    ! the circulation.
    !
    ! INPUTS
    !
    ! cl0: real(n_vort) -> cl0 of each 2D section
    ! cla: real(n_vort) -> lift curve slope of each 2D section [1/rad]

    implicit none

    ! Input variables
    integer, intent(in) :: n_vort
    real, intent(in) :: rho, cl0(n_vort), cla(n_vort), areas(n_vort)
    real, intent(in) :: Vlocal(3,n_vort), alphaLocal(n_vort), si(3,n_vort)

    ! Output variables
    real, intent(out) :: Fairf(3,n_vort)

    ! Working variables
    real :: Fairfm, cl, li(3), Uli(3), Vlocalm
    integer :: ii

    ! EXECUTION

    ! Loop over all horseshoe vortices
    do ii=1,n_vort

       ! Compute the magnitude of the force given by the 2D section
       ! L = 0.5*Vlocal²*area*cl
       cl = cl0(ii) + alphaLocal(ii)*cla(ii)
       call norm(Vlocal(:,ii), Vlocalm)
       Fairfm = 0.5*rho*Vlocalm*Vlocalm*areas(ii)*cl

       ! Get normalized lift direction
       call cross(Vlocal(:,ii),si(:,ii), li)
       call normalize(li, Uli)

       ! Store the force
       Fairf(:,ii) = Fairfm*Uli

    end do

  end subroutine compute_airfoil_forces

  !============================================================

  subroutine get_residuals(n_vort, X, Gama, twist, chords, cl0, cla, Vinf, rho, res)

    ! This subroutine executes all steps of the LLT method to compute
    ! residuals for a given Gama.

    implicit none

    ! Input variables
    integer, intent(in) :: n_vort
    real, intent(in) :: X(3,n_vort+1)
    real, intent(in) :: Gama(n_vort), twist(n_vort)
    real, intent(in) :: chords(n_vort), cl0(n_vort)
    real, intent(in) :: cla(n_vort), Vinf(3), rho

    ! Output variables
    real, intent(out) :: res(n_vort)

    ! Working variables
    integer :: ii
    real :: Uinf(3), AIC(3,n_vort,n_vort)
    real :: Uai(3,n_vort), Uni(3,n_vort), si(3,n_vort)
    real :: areas(n_vort), Vlocal(3,n_vort), alphaLocal(n_vort)
    real :: Fcirc(3,n_vort), Fairf(3,n_vort), deltaF(3)

    ! EXECUTION

    ! Normalize free-stream velocity
    call normalize(Vinf, Uinf)

    ! Get AIC matrix
    call get_AIC_matrix(n_vort, X, Uinf, AIC)

    ! Get geometric vectors
    call get_geom_vectors(n_vort, X, twist, chords, Uai, Uni, si, areas)

    ! Get induced velocities
    call get_local_vels(n_vort, Gama, AIC, Vinf, Uai, Uni, Vlocal, alphaLocal)

    ! Get circulation forces
    call compute_circulation_forces(n_vort, rho, Vlocal, si, Gama, Fcirc)

    ! Get airfoil forces
    call compute_airfoil_forces(n_vort, rho, cl0, cla, areas, Vlocal, alphaLocal, si, Fairf)

    ! The residual is the difference between the forces for each horseshoe vortex
    do ii=1,n_vort

       ! Get the difference in forces
       deltaF = Fcirc(:,ii) - Fairf(:,ii)

       ! The residual is the magnitude of the difference
       res(ii) = deltaF(1)**2 + deltaF(2)**2 + deltaF(3)**2

    end do

  end subroutine get_residuals

  !============================================================

  subroutine get_functions(n_vort, X, Gama, twist, chords, cl0, cla, Vinf, rho, Sref, CL, CD, L, D, Loads)

    ! This subroutine gives the functions of interest computed by the LLT.
    !
    ! INPUTS
    !
    ! n_vort: integer -> Number of horseshoe vortices.
    ! X: real(3,n_vort+1) -> Coordinates of the kinks of the horseshoe vortices.
    ! Gama: real(n_vort) -> Circulation intensity of each horseshoe vortex.
    ! twist: real(n_vort) -> Array of local incidence angles.
    ! chords: real(n_vort) -> Array of local chords.
    ! cl0: real(n_vort) -> cl0 of each 2D section
    ! cla: real(n_vort) -> lift curve slope of each 2D section [1/rad]
    ! Vinf: real(3) -> Free-stream velocity vector
    ! rho: real -> Air density
    !
    ! OUTPUTS
    !
    ! Sref: real -> Area of all panels
    ! CL: real -> Lift coefficient (adimensionalized by Sref)
    ! CD: real -> Drag coefficient (adimensionalized by Sref)
    ! L: real -> Lift force
    ! D: real -> Drag force
    ! Loads: real(n_vort) -> Total bending force over each bound vortex


    implicit none

    ! Input variables
    integer, intent(in) :: n_vort
    real, intent(in) :: X(3,n_vort+1)
    real, intent(in) :: Gama(n_vort), twist(n_vort)
    real, intent(in) :: chords(n_vort), cl0(n_vort)
    real, intent(in) :: cla(n_vort), Vinf(3), rho

    ! Output variables
    real, intent(out) :: Sref, CL, CD, L, D, Loads(n_vort)

    ! Working variables
    real :: Uinf(3), AIC(3,n_vort,n_vort)
    real :: Uai(3,n_vort), Uni(3,n_vort), si(3,n_vort)
    real :: areas(n_vort), Vlocal(3,n_vort), alphaLocal(n_vort)
    real :: Fcirc(3,n_vort), Fbody(3), Fbodym, Vinfm
    integer :: ii

    ! EXECUTION

    ! Normalize free-stream velocity
    call normalize(Vinf, Uinf)

    ! Get AIC matrix
    call get_AIC_matrix(n_vort, X, Uinf, AIC)

    ! Get geometric vectors
    call get_geom_vectors(n_vort, X, twist, chords, Uai, Uni, si, areas)

    ! Get induced velocities
    call get_local_vels(n_vort, Gama, AIC, Vinf, Uai, Uni, Vlocal, alphaLocal)

    ! Get forces on every bound vortex
    call compute_circulation_forces(n_vort, rho, Vlocal, si, Gama, Fcirc)

    ! Store bending forces
    Loads = Fcirc(3,:)

    ! Get total area and forces in the system of coordinates of the body
    Sref = 0.0
    Fbody = 0.0
    do ii = 1,n_vort
      Sref = Sref + areas(ii)
      Fbody(1) = Fbody(1) + Fcirc(1,ii)
      Fbody(2) = Fbody(2) + Fcirc(2,ii)
      Fbody(3) = Fbody(3) + Fcirc(3,ii)
    end do

    ! The drag is aligned with the free-stream
    call dot(Fbody, Uinf, D)

    ! The rest is the lift force
    call norm(Fbody, Fbodym)
    L = sqrt(Fbodym*Fbodym - D*D)

    ! Make results non dimensional
    call norm(Vinf, Vinfm)
    CL = 2.0*L/Sref/Vinfm/Vinfm/rho
    CD = 2.0*D/Sref/Vinfm/Vinfm/rho

  end subroutine get_functions

  !============================================================

  subroutine norm(a, am)

    implicit none

    real, intent(in) :: a(3)
    real, intent(out) :: am

    am = sqrt(a(1)*a(1) + a(2)*a(2) + a(3)*a(3))

  end subroutine norm

  !============================================================

  subroutine normalize(a, an)

    ! This subroutine normalizes a vector a to get the unitary vector an.

    implicit none

    real, intent(in) :: a(3)
    real, intent(out) :: an(3)
    real :: am

    call norm(a, am)
    an = a/am

  end subroutine normalize

  !============================================================

  subroutine dot(a, b, aDb)

    implicit none

    real, intent(in) :: a(3), b(3)
    real, intent(out) :: aDb

    aDb = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)

  end subroutine dot

  !============================================================

  subroutine cross(a, b, aCb)

    implicit none

    real, intent(in) :: a(3), b(3)
    real, intent(out) :: aCb(3)

    aCb(1) = a(2)*b(3) - a(3)*b(2)
    aCb(2) = a(3)*b(1) - a(1)*b(3)
    aCb(3) = a(1)*b(2) - a(2)*b(1)

  end subroutine cross

  !============================================================

  subroutine llt_main(n_vort, X, Gama, twist, chords, cl0, cla, Vinf, rho, res, Sref, CL, CD, L, D, Loads)

    ! Dummy subroutine to allow differentiation of get_residuals and
    ! get_functions in a single Tapenade call.
    !
    ! INPUTS
    !
    ! n_vort: integer -> Number of horseshoe vortices.
    ! X: real(3,n_vort+1) -> Coordinates of the kinks of the horseshoe vortices.
    ! Gama: real(n_vort) -> Circulation intensity of each horseshoe vortex.
    ! twist: real(n_vort) -> Array of local incidence angles.
    ! chords: real(n_vort) -> Array of local chords.
    ! cl0: real(n_vort) -> cl0 of each 2D section
    ! cla: real(n_vort) -> lift curve slope of each 2D section [1/rad]
    ! Vinf: real(3) -> Free-stream velocity vector
    ! rho: real -> Air density
    !
    ! OUTPUTS
    !
    ! res: real(n_vort) -> Residuals of the LLT method
    ! Sref: real -> Area of all panels
    ! CL: real -> Lift coefficient (adimensionalized by Sref)
    ! CD: real -> Drag coefficient (adimensionalized by Sref)
    ! L: real -> Total lift force
    ! D: real -> Total drag force
    ! Loads: real(n_vort) -> Total bending force over each bound vortex

    implicit none

    ! Input variables
    integer, intent(in) :: n_vort
    real, intent(in) :: X(3,n_vort+1)
    real, intent(in) :: Gama(n_vort), twist(n_vort)
    real, intent(in) :: chords(n_vort), cl0(n_vort)
    real, intent(in) :: cla(n_vort), Vinf(3), rho

    ! Output variables
    real, intent(out) :: res(n_vort), Sref, CL, CD, L, D
    real, intent(out) :: Loads(n_vort)

    ! EXECUTION

    ! Call relevant subroutines

    call get_residuals(n_vort, X, Gama, twist, chords, cl0, cla, Vinf, rho, res)
    call get_functions(n_vort, X, Gama, twist, chords, cl0, cla, Vinf, rho, Sref, CL, CD, L, D, Loads)

  end subroutine llt_main

end module llt_module
