subroutine Solver_SA          !the Solver
    implicit none
    integer::iter          !iterative variable
    integer::flag          !the variable to judge wheathe the mean density converges
    character(len=50)::filename
    
    write(*,*)  "Solver"
    
    call Grid
    call Allocate_memory  
    call Allocate_memory_SA
    call Flow_init_SA
    
    call Reorder   !reordering is needed for LU-SGS
    
    do iter = 1,itermax
        write(*,*)  iter
        
        call Mean_edge       !calculate the primative variables on every edge
        !turbulence
        muL =   1.45*( U(5,:)/R/U(1,:) )**(3.0/2)/ ( U(5,:)/R/U(1,:) + 110.0 ) * 1.0E-6
        muL_av = 1.45* U_av(6,:)**(3.0/2)/( U_av(6,:)+110.0 ) *1.0E-6
       
        call getMuT_SA( U,U_av,muL,muT,muT_av)
        
        if( schemeNum=="1" .OR. schemeNum=="3" ) then
            call Con_flux_JST        !calculate the convective flux using Roe scheme
            call Art_dissipation 
            dt = CFL*vol/(lamda_c + C_time*lamda_v )      !calculate the time step
        else
            call con_flux_Roe
            call Time_step
        end if
        
        call Gradient        !calculate the gradient of every cell using the Mean_edge subroutine
        
        call Resi_SA( U,U_av, alf,mul,mul_av ,Grad )      
        
        call Vis_flux        !calculate the viscous flux
        
        Rsi = Fc - Dissi - Fv
        
        call LU_SGS_SA
        
        call Converge(flag)
        
        if(flag == 1)  then
            write(*,*)  "...Converge..."
            exit 
       
        end if    
    end do
    
    filename ="steady"
    call Output_SA(filename)
    !call Cp_cal(filename)
    
end subroutine 
