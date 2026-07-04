module fem_module

  ! Module that holds necessary subroutines for
  ! the finite element method using beams.
  ! ATTENTION: This code assumes that the beam is horizontal and
  ! is only subject to vertical (bending) loads.

  implicit none

  real :: zero = 0.0
  real :: pi = 3.1415926535897932384626
  real :: eps = epsilon(zero)

contains

  subroutine beam_stiffness_matrix(X1i, X2i, ri, ti, E, Ki, Li, Ii)

    ! This subroutine creates the 4x4 matrix of a circular beam
    ! with 2 DOF per node.
    !
    ! INPUTS
    !
    ! X1i: real(3) -> Coordinates of the first node of the beam
    ! X2i: real(3) -> Coordinates of the second node of the beam
    ! ri: real -> Radius of the beam
    ! ti: real -> Wall thickness of the beam
    ! E: real -> Young's module of the beam material
    !
    ! OUTPUTS
    !
    ! Ki: real(4,4) -> Stiffness matrix of the beam
    ! Li: real -> Length of the beam
    ! Ii: real -> Moment of inertia of the beam

    implicit none

    ! Input variables
    real, intent(in) :: X1i(3), X2i(3), ri, ti, E

    ! Output variables
    real, intent(out) :: Ki(4,4), Li, Ii

    ! Working variables
    real :: kk

    ! EXECUTION

    ! Compute the length of the beam
    Li = X2i(2) - X1i(2)

    ! Compute the moment of inertia of a circular beam
    Ii = pi*ri**3*ti

    ! Compute the bending term that is common to all elements
    ! of the stiffness matrix
    kk = E*Ii/Li**3

    ! Populate the stiffness matrix
    Ki(1,1) = 12*kk
    Ki(1,2) = 6*Li*kk
    Ki(1,3) = -12*kk
    Ki(1,4) = 6*Li*kk

    Ki(2,1) = 6*Li*kk
    Ki(2,2) = 4*Li**2*kk
    Ki(2,3) = -6*Li*kk
    Ki(2,4) = 2*Li**2*kk

    Ki(3,1) = -12*kk
    Ki(3,2) = -6*Li*kk
    Ki(3,3) = 12*kk
    Ki(3,4) = -6*Li*kk

    Ki(4,1) = 6*Li*kk
    Ki(4,2) = 2*Li**2*kk
    Ki(4,3) = -6*Li*kk
    Ki(4,4) = 4*Li**2*kk

  end subroutine beam_stiffness_matrix

  !============================================================

  subroutine get_K_matrix(n_beams, X, r, t, E, K, n_dofs)

    ! This subroutine computes the global stiffness matrix of the problem.
    ! X is the list of beam nodes. We assume that the nodes in X are
    ! sequentially connected, that is, the first beam goes from X(:,1) to
    ! X(:,2), the second bar goes from X(:,2) to X(:,3), and so on.
    !
    ! INPUTS
    !
    ! n_beams: integer -> Number of beams.
    ! X: real(3,n_beams+1) -> Coordinates of the beam nodes.
    ! r: real(n_beams) -> Radius of each beam element.
    ! t: real(n_beams) -> Wall thickness of each beam element.
    ! E: real -> Young's module of the beam material
    !
    ! OUTPUTS
    !
    ! K: real(2*(n_beams+1),2*(n_beams+1)) -> Full stiffness matrix, including
    !                                         the constrained nodes.
    ! n_dofs: integer -> number of degree of freedom (size of Kf)

    implicit none

    ! Input variables
    integer, intent(in) :: n_beams
    real, intent(in) :: X(3,n_beams+1)
    real, intent(in) :: r(n_beams), t(n_beams), E

    ! Output variables
    real, intent(out) :: K(2*(n_beams+1),2*(n_beams+1))
    integer, intent(out) :: n_dofs

    ! Working variables
    integer :: jj
    real :: X1i(3), X2i(3), Ki(4,4), Li, Ii

    ! EXECUTION

    ! Initialize the stifness matrix
    K = zero

    ! Loop over all beams
    do jj=1,n_beams

       ! Get nodal coordinates of the current beam
       X1i = X(:,jj)
       X2i = X(:,jj+1)

       ! Compute the element stiffness matrix
       call beam_stiffness_matrix(X1i, X2i, r(jj), t(jj), E, Ki, Li, Ii)

       ! Assign the element matrix to the global matrix
       K(2*jj-1:2*jj+2,2*jj-1:2*jj+2) = K(2*jj-1:2*jj+2,2*jj-1:2*jj+2) + Ki

    end do

    ! Compute number of degrees of freedom
    n_dofs = 2*(n_beams+1)

  end subroutine get_K_matrix

  !============================================================

  subroutine apply_cons(n_dofs, n_cons, conIDs, K, f)

    ! This subroutine modifies the rows and columns of K and f to enforce
    ! the constraints from conIDs. For instance, if we want to set the i-th DOF
    ! to zero, we set f(i)=0, K(i,:)=0, K(:,i)=1, and K(i,i)=1.
    !
    ! INPUTS
    !
    ! n_dofs: integer -> Number of degrees of freedom.
    ! n_cons: integer -> Number of constrained DOFs
    ! conIDs: integer(n_cons) -> Indices of the DOFs of the Kf matrix that
    !                            should be constrained.
    !
    ! OUTPUTS
    !
    ! INPUTS-OUTPUTS
    !
    ! K: real(n_dofs,n_dofs) -> Full stiffness matrix, including
    !                           the constrained nodes.
    ! f: real(n_dofs) -> Full force vector, including the
    !                    constrained nodes.

    implicit none

    ! Input variables
    integer, intent(in) :: n_dofs, n_cons
    integer, intent(in) :: conIDs(n_cons)

    ! Output variables
    real, intent(inout) :: K(n_dofs,n_dofs)
    real, intent(inout) :: f(n_dofs)

    ! Working variables
    integer :: ii, curr_con

    ! EXECUTION

    ! First we will create the array freeIDs containing all indices
    ! of the free DOFs. Then we will use this array to crop Kf and ff.

    ! Now we set the indices of the constrained DOFs to zero
    do ii=1,n_cons
       curr_con = conIDs(ii)

       ! We will modify K and f so that d(curr_con) is forced to be zero.
       ! We can do that by setting:
       ! f(curr_con)=0, K(curr_con,:)=0, K(:,curr_con)=1,
       ! and K(curr_con,curr_con)=1.

       ! Modify the K matrix
       K(curr_con,:) = 0.0
       K(:,curr_con) = 0.0
       K(curr_con,curr_con) = 1.0

       ! Modify the f vector
       f(curr_con) = 0.0

    end do

  end subroutine apply_cons

  !============================================================

  subroutine get_residuals(n_beams, n_cons, X, r, t, f, d, E, conIDs, res)

    ! This subroutine computes the residuals of the FEM problem.
    !
    ! INPUTS
    !
    ! n_beams: integer -> Number of beams.
    ! n_cons: integer -> Number of constrained DOFs
    ! X: real(3,n_beams+1) -> Coordinates of the beam nodes.
    ! r: real(n_beams) -> Radius of each beam element.
    ! t: real(n_beams) -> Wall thickness of each beam element.
    ! f: real(2*(n_beams+1)) -> Forces on all nodes (including
    !                           the constrained ones).
    ! d: real(2*(n_beams+1)) -> Displacements values of all DOFs (including
    !                           the constrained ones).
    ! E: real -> Young's module of the beam material
    ! conIDs: integer(n_cons) -> Indices of the DOFs of the Kf matrix that
    !                            should be constrained.
    !
    ! OUTPUTS
    !
    ! res: real(2*(n_beams+1)) -> Residuals of the FEM model
    !                             considering all DOFs (including the
    !                             constrained ones). The residuals
    !                             of the constrained DOFs are the
    !                             displacements themselves, so that the
    !                             solver drives them to zero. We
    !                             modify the K matrix in appy_cons to
    !                             achieve that.

    implicit none

    ! Input variables
    integer, intent(in) :: n_beams, n_cons
    real, intent(in) :: X(3,n_beams+1)
    real, intent(in) :: r(n_beams), t(n_beams)
    real, intent(in) :: f(2*(n_beams+1)), d(2*(n_beams+1)), E
    integer, intent(in) :: conIDs(n_cons)

    ! Output variables
    real, intent(out) :: res(2*(n_beams+1))

    ! Working variables
    real :: K(2*(n_beams+1), 2*(n_beams+1)), fcopy(2*(n_beams+1))
    real :: Kd(2*(n_beams+1))
    integer :: n_dofs

    ! EXECUTION

    ! Build the full stiffness matrix
    call get_K_matrix(n_beams, X, r, t, E, K, n_dofs)

    ! Make a copy of the force vector so we can modify its elements
    fcopy = f

    ! Apply the constraints
    call apply_cons(n_dofs, n_cons, conIDs, K, fcopy)

    ! Compute K*d
    call matdotvec(n_dofs, K, d, Kd)

    ! Compute the residuals
    res = Kd - fcopy

  end subroutine get_residuals

  !============================================================

  subroutine get_functions(n_beams, n_cons, X, r, t, f, d, E, rho, sigmaY, &
                           pKS, conIDs, mass, margins, KSmargin)

    ! This subroutine computes the residuals of the FEM problem.
    !
    ! INPUTS
    !
    ! n_beams: integer -> Number of beams.
    ! n_cons: integer -> Number of constrained DOFs
    ! X: real(3,n_beams+1) -> Coordinates of the beam nodes.
    ! r: real(n_beams) -> Radius of each beam element.
    ! t: real(n_beams) -> Wall thickness of each beam element.
    ! f: real(2*(n_beams+1)) -> Forces on all nodes (including
    !                           the constrained ones).
    ! d: real(2*(n_beams+1)) -> Displacements values of all DOFs (including
    !                           the constrained ones).
    ! E: real -> Young's module of the beam material.
    ! rho: real -> Density of the beam material.
    ! sigmaY: real -> Yield stress of the beam material.
    ! pKS: real -> Constant of the KS functons used to aggregate stresses.
    ! conIDs: integer(n_cons) -> Indices of the DOFs of the Kf matrix that
    !                            should be constrained.
    !
    ! OUTPUTS
    !
    ! mass: real -> Mass of the entire structure.
    ! margins: real(2*n_beams) -> Margin from Von Mises stress at each node
    !                             with respect to the yield stress.
    !                             The margins are computed as:
    !                             margin = 1.0 - sigma_vm/sigmaY
    !                             Negative margins indicate failure.
    !                             Internal nodes have two values
    !                             corresponding to the beam sections on
    !                             either side. KSmargin is the aggregation of
    !                             these values.
    ! KSmargin: real -> Aggregated safety margin. This value
    !                   is a conservative estimate of the minimum
    !                   margin value in the structure. This estimate
    !                   gets closer to the real value as we increase
    !                   pKS, but this may lead to numerical issues.

    implicit none

    ! Input variables
    integer, intent(in) :: n_beams, n_cons
    real, intent(in) :: X(3,n_beams+1)
    real, intent(in) :: r(n_beams), t(n_beams)
    real, intent(in) :: f(2*(n_beams+1)), d(2*(n_beams+1))
    real, intent(in) :: E, rho, sigmaY, pKS
    integer, intent(in) :: conIDs(n_cons)


    ! Output variables
    real, intent(out) :: mass, margins(2*n_beams), KSmargin

    ! Working variables
    real :: X1i(3), X2i(3), Ki(4,4), Li, Ii, di(4), fi(4), Ai
    real :: sigma1, sigma2, tau, sigma_vm1, sigma_vm2
    real :: margin_min, margin_sum
    integer :: jj

    ! EXECUTION

    ! Initialize the mass value
    mass = 0

    ! Loop over every element to compute contributions of mass and stress
    do jj=1,n_beams

       ! Get nodal coordinates of the current beam
       X1i = X(:,jj)
       X2i = X(:,jj+1)

       ! Get element properties
       call beam_stiffness_matrix(X1i, X2i, r(jj), t(jj), E, Ki, Li, Ii)

       ! Compute section area
       Ai = 2.0*pi*r(jj)*t(jj)

       ! Compute mass contribution
       mass = mass + rho*Ai*Li

       ! Get DOF values of the current nodes
       di = d(2*jj-1:2*jj+2)

       ! Compute corresponding local nodal forces
       call matdotvec(4, Ki, di, fi)

       ! Compute bending stress on node 1 (sigma = M*y/I)
       sigma1 = fi(2)*r(jj)/Ii

       ! Compute bending stress on node 2 (sigma = M*y/I)
       sigma2 = fi(4)*r(jj)/Ii

       ! Compute shear stress on the beam (it is constant)
       ! We need 2.0 because shear is maximum at the midplane
       tau = 2.0*fi(1)/Ai

       ! Compute the Von Mises stress for each node
       ! This is a conservative measure, since we are using the
       ! maximum bending stress and the maximum shear stress, even though
       ! they are at different locations of the cross-section.
       sigma_vm1 = sqrt(sigma1**2 + 3.0*tau**2)
       sigma_vm2 = sqrt(sigma2**2 + 3.0*tau**2)

       ! Compute the safety margin with respect to the failure
       margins(2*jj-1) = 1.0 - sigma_vm1/sigmaY
       margins(2*jj) = 1.0 - sigma_vm2/sigmaY

    end do

    ! The KS function estimates the maximum value from the set.
    ! But the minimum margin is the most critical one.
    ! So we flip the sign of the margins to
    ! compute the KS value.

    ! Find the minimum safety margin (most critical case)
    ! NOTE: I deactivated this since the changes in the minimum values
    ! caused discrepancies with respect to the finite difference test
    ! when checking derivatives.
    margin_min = 0.0 !1.0
    !do jj=1,2*n_beams
    !   margin_min = min(margin_min, margins(jj))
    !end do

    ! Compute sum of the KS expression inside the natural logarithm.
    ! Note that we flipped the margins signs.
    margin_sum = 0.0
    do jj=1,2*n_beams
      margin_sum = margin_sum + exp(pKS*(-margins(jj) + margin_min))
    end do

    ! Aggregate stresses
    KSmargin = margin_min - 1.0/pKS*log(margin_sum)

  end subroutine get_functions

  !============================================================

  subroutine matdotvec(n, K, d, Kd)
    ! This subroutine does the dot product between matrix K
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

  !============================================================

  subroutine fem_main(n_beams, n_cons, X, r, t, f, d, &
                      E, rho, sigmaY, pKS, conIDs, &
                      res, mass, margins, KSmargin)

    ! Dummy subroutine to allow differentiation of get_residuals and
    ! get_functions in a single Tapenade call.
    !
    ! INPUTS
    !
    ! n_beams: integer -> Number of beams.
    ! n_cons: integer -> Number of constrained DOFs
    ! X: real(3,n_beams+1) -> Coordinates of the beam nodes.
    ! r: real(n_beams) -> Radius of each beam element.
    ! t: real(n_beams) -> Wall thickness of each beam element.
    ! f: real(2*(n_beams+1)) -> Forces on all nodes (including
    !                           the constrained ones).
    ! d: real(2*(n_beams+1)) -> Displacements values of all DOFs (including
    !                           the constrained ones).
    ! E: real -> Young's module of the beam material.
    ! rho: real -> Density of the beam material.
    ! sigmaY: real -> Yield stress of the beam material.
    ! pKS: real -> Constant of the KS functons used to aggregate stresses.
    ! conIDs: integer(n_cons) -> Indices of the DOFs of the Kf matrix that
    !                            should be constrained.
    !
    ! OUTPUTS
    !
    ! res: real(2*(n_beams+1)) -> Residuals of the FEM model
    !                             considering all DOFs (including the
    !                             constrained ones). The residuals
    !                             of the constrained DOFs are the
    !                             displacements themselves, so that the
    !                             solver drives them to zero. We
    !                             modify the K matrix in appy_cons to
    !                             achieve that.
    ! mass: real -> Mass of the entire structure.
    ! margins: real(2*n_beams) -> Margin from Von Mises stress at each node
    !                             with respect to the yield stress.
    !                             The margins are computed as:
    !                             margin = 1.0 - sigma_vm/sigmaY
    !                             Negative margins indicate failure.
    !                             Internal nodes have two values
    !                             corresponding to the beam sections on
    !                             either side. KSmargin is the aggregation of
    !                             these values.
    ! KSmargin: real -> Aggregated safety margin. This value
    !                   is a conservative estimate of the minimum
    !                   margin value in the structure. This estimate
    !                   gets closer to the real value as we increase
    !                   pKS, but this may lead to numerical issues.

    implicit none

    ! Input variables
    integer, intent(in) :: n_beams, n_cons
    real, intent(in) :: X(3,n_beams+1)
    real, intent(in) :: r(n_beams), t(n_beams)
    real, intent(in) :: f(2*(n_beams+1)), d(2*(n_beams+1))
    real, intent(in) :: E, rho, sigmaY, pKS
    integer, intent(in) :: conIDs(n_cons)


    ! Output variables
    real, intent(out) :: res(2*(n_beams+1)), mass, margins(2*n_beams), KSmargin

    ! EXECUTION

    ! Call relevant subroutines

    call get_residuals(n_beams, n_cons, X, r, t, f, d, E, conIDs, res)

    call get_functions(n_beams, n_cons, X, r, t, f, d, &
                       E, rho, sigmaY, pKS, conIDs, &
                       mass, margins, KSmargin)

  end subroutine fem_main

end module fem_module
