module mymodule

contais

    subroutine average(nvec, vec, avg)
    
    	implicit none
    	
    	! Input variables
    	integer, intent(in) :: nvec
    	real, intent(in) :: vec(nvec)
    	
    	!Output variables
    	real, intent(out) :: avg
    	
    	!Working variables
    	real :: vectotal
    	
    	
    	
    	!EXECTION
    	
    	vectotal = sum(vec)
    	
    	avg = vectotal/nvec
    
    
    
    
    
    end subroutine average







end module mymodule
