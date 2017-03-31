subroutine flowInit
    implicit none
    
    if( moveNum== 0 )  then
        call flowInit_static  
    else
        call flowInit_dynamic
    end if
    
end subroutine
    