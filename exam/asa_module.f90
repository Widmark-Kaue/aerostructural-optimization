module asa_module

  ! Module that performs aerostructural analysis.
  ! This module also holds necessary subroutines for
  ! transfering loads and displacement betweenn LLT and FEM models.
  ! ATTENTION: This code assumes that the beam
  ! is only subject to vertical (bending) loads.

  implicit none

  real :: zero = 0.0
  real :: pi = 3.1415926535897932384626
  real :: eps = epsilon(zero)

contains

  subroutine transfer_disps(n_nodes, Xfem, d, Xllt)

    ! This subroutine takes the FEM nodes and displacements to generate
    ! the LLT nodes.
    !
    ! INPUTS
    !
    ! n_nodes: integer -> Number of nodes (dimension 2 of Xfem)
    ! Xfem: real(3,n_nodes) -> Coordinates of the FEM nodes
    ! d: real(2*n_nodes) -> Displacements and rotations of the FEM nodes,
    !                       according to the FEM formulation:
    !                       d = [d1, r1, d2, r2, ..., dn, rn]
    !                       where di are vertical displacements (z-dir)
    !                       and ri are rotations along the bending axis (x-dir)
    !
    ! OUTPUTS
    !
    ! Xllt: real(3,n_nodes) -> Coordinates of the deformed wing nodes, which
    !                          are used as kinks of the bound vortices.

    implicit none

    ! Input variables
    integer, intent(in) :: n_nodes
    real, intent(in) :: Xfem(3,n_nodes), d(2*n_nodes)

    ! Output variables
    real, intent(out) :: Xllt(3,n_nodes)

    ! EXECUTION

    ! Copy Xfem to Xllt, but just apply displacements in the z direction
    Xllt(1,:) = Xfem(1,:)
    Xllt(2,:) = Xfem(2,:)
    Xllt(3,:) = Xfem(3,:) + d(::2)

  end subroutine transfer_disps

  !============================================================

  subroutine transfer_forces(n_nodes, Xfem, Loads, f)

    ! This subroutine takes the vertical loads of the LLT elements and
    ! then computes the consistent forces at the FEM nodes.
    ! Here we use consistent forces for a constant distributed load
    ! over the LLT element.
    !
    ! INPUTS
    !
    ! n_nodes: integer -> Number of nodes (dimension 2 of Xfem)
    ! Xfem: real(3,n_nodes) -> Coordinates of the FEM nodes
    ! Loads: real(n_nodes-1) -> Total bending force over each bound vortex
    !
    ! OUTPUTS
    !
    ! f: real(2*n_nodes) -> Consistent forces at each FEM node,
    !                       according to the FEM formulation:
    !                       f = [f1, m1, f2, m2, ..., fn, mn]
    !                       where fi are vertical forces (z-dir)
    !                       and mi are bending moments (x-dir)

    implicit none

    ! Input variables
    integer, intent(in) :: n_nodes
    real, intent(in) :: Xfem(3,n_nodes), Loads(n_nodes-1)

    ! Output variables
    real, intent(out) :: f(2*n_nodes)

    ! Working variables
    real :: L(n_nodes-1)
    integer :: ii

    ! EXECUTION

    ! Compute the panel lengths along the span axis (y-dim)
    L = Xfem(2,2:n_nodes) - Xfem(2,1:n_nodes-1)

    ! Initialize consistent force vector
    f = zero

    ! Compute contributions of each horseshoe vortex to the consistent forces
    do ii=1,n_nodes-1
       f(2*ii-1) = f(2*ii-1) + Loads(ii)/2
       f(2*ii) = f(2*ii) + Loads(ii)*L(ii)/12
       f(2*ii+1) = f(2*ii+1) + Loads(ii)/2
       f(2*ii+2) = f(2*ii+2) - Loads(ii)*L(ii)/12
    end do

  end subroutine transfer_forces

  !============================================================

  subroutine asa_main(n_panels, Gama, twist, chords, &
                      span, r, t, d, &
                      cl0, cla, Vinf, rhoAir, &
                      E, rhoMat, sigmaY, pKS, &
                      CD0, fixedMass, g, Endurance, &
                      TSFC, loadFactor, &
                      resllt, resfem, liftExcess, &
                      margins, KSmargin, FB, Weight, Sref, CL)

    ! This subroutine computes the residuals of the coupled aerostructural
    ! problem.
    !
    ! INPUTS
    !
    ! n_panels: integer -> Number of horseshoe panels and beams (the same).
    ! Gama: real(n_panels) -> Circulation intensity of each horseshoe vortex.
    ! twist: real(n_panels) -> Array of local incidence angles.
    ! chords: real(n_panels) -> Array of local chords.
    ! span: real -> Wing span.
    ! r: real(n_panels) -> Radius of each beam element.
    ! t: real(n_panels) -> Wall thickness of each beam element.
    ! d: real(2*(n_panels+1)) -> Displacements values of all DOFs (including
    !                           the constrained ones).
    ! cl0: real(n_panels) -> cl0 of each 2D section
    ! cla: real(n_panels) -> lift curve slope of each 2D section [1/rad]
    ! Vinf: real(3) -> Free-stream velocity vector
    ! rhoAir: real -> Air density
    ! E: real -> Young's module of the beam material.
    ! rhoMat: real -> Density of the beam material.
    ! sigmaY: real -> Yield stress of the beam material.
    ! pKS: real -> Constant of the KS functons used to aggregate stresses.
    ! CD0: real -> Parasite drag that will be added to the induced drag.
    ! fixedMass: real -> Fixed mass of the aircraft less the wing structure.
    ! g: real -> Acceleration of gravity.
    ! Endurance: real -> Desired endurance of the aircraft.
    ! TSFC: real -> Thrust specific fuel consumption of the engine.
    ! loadFactor: real -> Factor applied to the structural failure margin
    !                     to consider a multiplicative factor in the loads.
    !
    ! OUTPUTS
    !
    ! resllt: real(n_panels) -> Residuals of the LLT method
    ! resfem: real(2*(n_panels+1)) -> Residuals of the FEM model
    !                                 considering all DOFs (including the
    !                                 constrained ones). The residuals
    !                                 of the constrained DOFs are the
    !                                 displacements themselves, so that the
    !                                 solver drives them to zero. We
    !                                 modify the K matrix in appy_cons to
    !                                 achieve that.
    ! liftExcess: real -> Difference between the generated lift and the.
    !                     aircraft weight. This value is normalized:
    !                     liftExcess = Lift/Weigth - 1
    ! margins: real(2*n_panels) -> Margin from Von Mises stress at each node
    !                              with respect to the yield stress.
    !                              The margins are computed as:
    !                              margin = 1.0 - sigma_vm/sigmaY
    !                              Negative margins indicate failure.
    !                              Internal nodes have two values
    !                              corresponding to the beam sections on
    !                              either side. KSmargin is the aggregation of
    !                              these values.
    ! KSmargin: real -> Aggregated safety margin. This value
    !                   is a conservative estimate of the minimum
    !                   margin value in the structure. This estimate
    !                   gets closer to the real value as we increase
    !                   pKS, but this may lead to numerical issues.
    ! FB: real -> Fuel burn (in Force units).
    ! Weight: real -> Total weight of the aircraft (in Force units)
    !                 (fixed weight + structure weight + fuel weight)
    ! Sref: real -> Area of all panels
    ! CL: real -> Lift coefficient (adimensionalized by Sref)

    use llt_module, only : llt_main
    use fem_module, only : fem_main
    implicit none

    ! Input variables
    integer, intent(in) :: n_panels
    real, intent(in) :: Gama(n_panels), twist(n_panels)
    real, intent(in) :: chords(n_panels)
    real, intent(in) :: span
    real, intent(in) :: r(n_panels), t(n_panels)
    real, intent(in) :: d(2*(n_panels+1))
    real, intent(in) :: cl0(n_panels), cla(n_panels), Vinf(3), rhoAir
    real, intent(in) :: E, rhoMat, sigmaY, pKS
    real, intent(in) :: CD0, fixedMass, g, Endurance, TSFC, loadFactor

    ! Output variables
    real, intent(out) :: resllt(n_panels), resfem(2*(n_panels+1))
    real, intent(out) :: liftExcess, margins(2*n_panels)
    real, intent(out) :: KSmargin, FB, Weight, Sref, CL

    ! Working variables
    real :: Xfem(3,n_panels+1)
    real :: Xllt(3,n_panels+1), CD, Lift, Drag
    real :: Loads(n_panels), f(2*(n_panels+1)), structMass
    real :: Vinfm2
    integer :: n_cons, conIDs(2)
    integer :: ii, n_nodes

    ! EXECUTION

    ! Identify the indices of the structural DOFs that should
    ! be constrained, which corresponds to the central node.
    ! The first node of the wing has degrees of freedom 1 and 2,
    ! the second node has degrees of freedom 3 and 4, and so on.
    ! If we have n_panels elements, the degrees of freedom of the
    ! central node will be n_panels+1 and n_panels+2.
    n_cons = 2
    conIDs(1) = n_panels+1
    conIDs(2) = n_panels+2
    n_nodes = n_panels + 1

    ! ==================
    ! ADD YOUR CODE HERE
    ! ==================


  end subroutine asa_main

end module asa_module
