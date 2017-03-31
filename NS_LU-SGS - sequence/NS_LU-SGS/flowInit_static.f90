subroutine flowInit_static
    implicit none
    
    select case( turModelNum )
        case( 0 ) 
            call flowInit_turNone_static
        case( 1 )
            call flowInit_SA_static
        case( 2 )
            call flowInit_SST_static
        case default
      
    end select
    
end subroutine
    
subroutine flowInit_turNone_static

    implicit none
   
    U(1,:)=rou_inf
    U(2,:)=u_inf
    U(3,:)=v_inf
    U(5,:)=p_inf
    
    W(1,:)=rou_inf
    W(2,:)=rou_inf*u_inf
    W(3,:)=rou_inf*v_inf
    W(5,:)=p_inf/(gamma-1.0) + rou_inf*(u_inf**2+v_inf**2)/2.0 
    
    muT =0.0
    muT_av = 0.0

end subroutine
    
    
subroutine flowInit_SA_static      !initialize the flow field

    implicit none
    
    U(1,:)=rou_inf
    U(2,:)=u_inf
    U(3,:)=v_inf
    U(5,:)=p_inf
    
    W(1,:)=rou_inf
    W(2,:)=rou_inf*u_inf
    W(3,:)=rou_inf*v_inf
    W(5,:)=p_inf/(gamma-1.0) + rou_inf*(u_inf**2+v_inf**2)/2.0 
   
    !initilize the turbulent flow
    U_SA = 0.1*muL_inf/rou_inf
    W_SA = rou_inf*U_SA
    
    muT = 0.1*muL_inf
   
end subroutine

subroutine flowInit_SST_static      !initialize the flow field

    implicit none

    U(1,:)=rou_inf
    U(2,:)=u_inf
    U(3,:)=v_inf
    U(5,:)=p_inf
    
    W(1,:)=rou_inf
    W(2,:)=rou_inf*u_inf
    W(3,:)=rou_inf*v_inf
    W(5,:)=p_inf/(gamma-1.0) + rou_inf*(u_inf**2+v_inf**2)/2.0    
    
    !initilize the turbulent flow

  
    U_SST(2,:) = C1_SST*sqrt(u_inf**2 + v_inf**2)/ L_SST
    U_SST(1,:) = muL_inf*10.0**(-C2_SST)/rou_inf *U_SST(2,:)
    
    W_SST(1,:) = U(1,:) * U_SST(1,:)
    W_SST(2,:) = U(1,:) * U_SST(2,:)
    
    muT = muL_inf*10.0**(-C2_SST)
    
end subroutine
