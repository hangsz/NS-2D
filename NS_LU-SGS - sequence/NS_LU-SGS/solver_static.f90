subroutine solver_static          !the Solver
    implicit none
    integer::iter          !iterative variable
    integer::flag          !the variable to judge wheathe the mean density converges
    character(len=50)::filename
    real(8)::oldRsi(5)
    !--------------------------------
    write(*,*)  "solver_static"
    
    do iter = 1,itermax
        write(*,*)  iter
        
        call meanEdge
        call gradient
        
        !turbulence
         !----------------------------------
        ! calculate the turbulent viscosity
        select case( turModelNum ) 
        case( 0 )     
        case( 1 )  
            call getMuT_SA(U_av,muT,muT_av,U_Rot)
            call Resi_SA( U,U_av,muL,muL_av ,Grad ,U_Rot )  
            call LU_SGS_SA   
        case( 2 )
            call getMUT_SST(U_av,muT,muT_av,U_Rot)
            call Resi_SST( u_inf,v_inf,U,U_av, muL,muL_av,muT,muT_av,Grad,U_Rot )     
            call LU_SGS_SST   
        case default 
        end select     
        !--------------------------------------------
        ! renew the turbulent viscosity on the edge
        select case( turModelNum ) 
        case( 1 )
            call getMuT_SA(U_av,muT,muT_av,U_Rot)
        case( 2 )
            call getMUT_SST(U_av,muT,muT_av,U_Rot)
        case default 
        end select
            
        select case( conFluxNum )
        case( 1 )
            call conFlux_JST        !calculate the convective flux using Roe scheme
            call artDissipation 
        case( 2 ) 
            Dissi = 0.0
            call conFlux_Roe
            call timeStep
        end select
        
        call visFlux        !calculate the viscous flux
        
        Rsi = Fc - Dissi - Fv
        
        call LU_SGS
        !-----------------------------------
        if( iter == 1 )  oldRsi = sum( abs(Rsi(1,:)))/ncells
                
        call Converge(flag,oldRsi)
        
        if(flag == 1)  then
            write(*,*)  "...Converge..."
            exit 
        end if  
           
    end do
    filename ="steady"
    call Output(filename)   
     
end subroutine 
