program main
    use llt, only: llt_main
    
    implicit none
    


    ! Input variables
    !INPUTS
    integer :: ii
    real :: twist(6) = (/ 0., 0., 0., 0., 0., 0. /)
    real :: gama(6) = (/ 1., 1., 1., 1., 1., 1. /)
    real :: span = 8.0
    real :: chord = 1.0
    real :: cl0 = 0.0
    real :: cla = 6.283185307179586
    real :: alpha = 0.08726646259971647
    real :: vinf = 10
    real :: rho_air = 1.225
    real :: res_llt(6)
    real :: CL, CD
    call llt_main(6, twist, gama, span, chord, cl0, cla, alpha, vinf, rho_air, res_llt, CL, CD)

    print *, "LLT residuals at each wing section: ["

    do ii = 1, 6
        print *, res_llt(ii)
    end do
    print *, "]"
    
    print *, "Total lift coefficient of the wing: ", CL
    print *, "Total drag coefficient of the wing: ", CD
end program  main