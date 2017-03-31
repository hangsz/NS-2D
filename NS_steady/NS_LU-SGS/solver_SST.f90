subroutine Solver_SST          !the Solver
    implicit none
    integer::i,j
    integer::iter          !iterative variable
    integer::count
    integer::flag          !the variable to judge wheathe the mean density converges
    character(len=50)::filename
    
    write(*,*)  "Solver_SST"
    
    call Grid
    call Allocate_memory  
    call Allocate_memory_SST
    
    call Flow_init_SST
    
    call Reorder   !reordering is needed for LU-SGS
       
    !-------------------------------------------
    !geometry set
    !translate the grid to where the rotational center's coordinates are zero
    
         
    do iter =1,itermax
        write(*,*)   iter
        
        call Mean_edge       !calculate the primative variables on every edge
              
        !-------------------------
        !the viscosity on the edges
        muL =   1.45*( U(5,:)/R/U(1,:) )**(3.0/2)/ ( U(5,:)/R/U(1,:) + 110.0 ) * 1.0E-6
        muL_av = 1.45*U_av(6,:)**(3.0/2)/( U_av(6,:)+110.0 ) *1.0E-6
        !turbulence
        
        call getMUT_SST(U,U_av,muL,muT,muT_av)
        !-----------------
        call Mean_edge_SST( u_inf,v_inf, U_av,muL_av)   ! its needed to calculate convective flux
        
        if( schemeNum=="2" .OR. schemeNum=="4" ) then
            call Con_flux_JST        ! calculate the convective flux using Roe scheme
            call Art_dissipation 
            dt = CFL*vol/(lamda_c + C_time*lamda_v )       !calculate the time step
        else
            call con_flux_Roe
            call time_step
        end if
        
        
        call Gradient            !calculate the gradient of every cell using the Mean_edge subroutine
           
        call Resi_SST( U,U_av, alf,muL,muL_av,muT,muT_av,Grad )     !calculate the muT
    
        call Vis_flux        !calculate the viscous flux
        
        Rsi = Fc - Dissi - Fv
                          
        call LU_SGS_SST
        
        call Converge(flag)
        
        if(flag == 1)  then
            write(*,*)  "...Converge..."
            exit 
        end if    
    end do
          
    filename ="steady"
    call Output_SST(filename)
    
    !write(*,*)  "Cp_cal"
    !call Cp_cal(filename)
        
end subroutine 