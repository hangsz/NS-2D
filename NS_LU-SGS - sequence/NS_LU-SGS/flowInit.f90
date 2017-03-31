subroutine flowInit
    implicit none
    
    a_inf=sqrt(p_inf*gamma/rou_inf)
    u_inf=Ma_inf*a_inf*cosd(att)        
    v_inf=Ma_inf*a_inf*sind(att)
   
    U_Rot = 0.0
    
    if( moveNum== 0 )  then
        call flowInit_static  
    else
        call flowInit_dynamic
    end if
    
end subroutine
    