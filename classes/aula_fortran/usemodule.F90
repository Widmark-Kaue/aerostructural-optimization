program main

    use mymodule, only : average
    
    implicit none
    
    integer :: nvec
    real, dimension(:), allocatable :: vec
    real :: avg
    
    nvec = 5
    
    allocate(vec(nvec))
    
    do ii = 1,nvec
    	vec(ii) = ii*1.5
    end do
    
    print*,'vec:',vec
    
    call average(nvec, vec, avg)
    
    print *,'avg',avg
    
    deallocate(vec)







end main
