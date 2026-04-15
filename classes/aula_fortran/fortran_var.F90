program main	

    implicit none
    
    integer :: ii, jj 
    real :: aa, bb
    real, dimension(4) :: vec1
    real ::vec2(4)
    real :: MM(3,2)
    real, dimension (:,:), allocatable :: PP
    character(len=15) :: string
    
    
    ! EXECUTION 
    
    ii = 123
    jj = 321
    
    print *, 'Integers'
    print *, ii
    print *, jj 
    
    aa = 3.0
    bb = 4.0
    print *, 'Reais'
    print *,sqrt(aa**2+bb**2)
    
    vec1 = (/ 0.5, 1.5, 2.5, 3.5/)
    print *, 'Vectors'
    print *, vec1
    print *, vec1(2)
    print *, vec1(2:4)
    
    MM = 2.0
    print *, 'Matrix'
    print *, MM
    MM(2,1) = 0.0
    print *, 'Matrix2'
    
    do jj = 1,3
    
       print*,MM(jj,:)
    
    end do
    
    aa = 2.5
    
    print *, 'If example'
    
    if (aa <=3.0) then
    	print *, 'aa menor que 3.0'
    else if (aa <= 4.0) then
    	print *, 'aa entre 3.0 e 4.0'
    else
    	print *, 'aa maior que 4.0'
    end if

    print *, size(vec1)
    ii = 5
    jj = 2
    allocate(PP(ii,jj))
    
    PP = 1.0
    
    print *, 'Allocated'
    do ii = 1,5
    	print *,PP(ii,:)
    end do
    
    deallocate(PP)
    
    string = 'test'
    print *, 'String'
    print *,string
    
end program main
