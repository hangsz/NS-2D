subroutine solver
    implicit none
    
    call Grid
    call allocateMemory  
    
    call Reorder   !reordering is needed for LU-SGS
    
    select case( turModelNum ) 
        case( 1 )
            call allocateMemory_SA
        case( 2 )
            call allocateMemory_SST
        case default 
    end select
             
    call flowInit
    
    if( moveNum== 0 )  then
        call solver_static  
    else
        call outputFreeFlow
    
        call solver_dynamic
    end if
    
    
end subroutine
    

    
    