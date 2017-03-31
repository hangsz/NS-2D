subroutine solver_static          !the Solver
    implicit none
    integer::iter          !iterative variable
    integer::flag          !the variable to judge wheathe the mean density converges
    character(len=50)::filename
    
    write(*,*)  "solver_static"
    
    
    
    do iter = 1,itermax
        write(*,*)  iter
        
        call meanEdge       !calculate the primative variables on every edge
        !turbulence
        muL =   1.45*( U(5,:)/R/U(1,:) )**(3.0/2)/ ( U(5,:)/R/U(1,:) + 110.0 ) * 1.0E-6
        muL_av = 1.45* U_av(6,:)**(3.0/2)/( U_av(6,:)+110.0 ) *1.0E-6
       
        select case( turModelNum ) 
            case( 1 )
                call getMuT_SA( U,U_av,muL,muT,muT_av,U_Rot)
                call meanEdge_SA( U_av,U_Rot)  
            case( 2 )
                
                call getMUT_SST(U,U_av,muL,muT,muT_av,U_Rot)
                call meanEdge_SST( u_inf,v_inf, U_av,muL_av, U_Rot)   ! its needed to calculate convective flux 
            case default 
            
        end select
            
        select case( conFluxNum )
        case( 1 )
            call conFlux_JST        !calculate the convective flux using Roe scheme
            call artDissipation 
            dt = CFL*vol/(lamda_c + C_time*lamda_v )      !calculate the time step
        case( 2 ) 
            call conFlux_Roe
            call timeStep
        end select
        

        call gradient        !calculate the gradient of every cell using the Mean_edge subroutine 
        
        call visFlux        !calculate the viscous flux
        
        Rsi = Fc - Dissi - Fv
        
        select case( turModelNum ) 
            case( 0 ) 
                call LU_SGS
            case( 1 )
                
                call Resi_SA( U,U_av, alf,muL,muL_av ,Grad ,U_Rot )  
                call LU_SGS_SA
                
            case( 2 )
                
                call Resi_SST( U,U_av, alf,muL,muL_av,muT,muT_av,Grad,U_Rot )     
                call LU_SGS_SST
                
            case default 
            end select       
        
            if( iter == 1 )  oldRsi = sum( abs(Rsi(1,:)))/ncells
                
            call Converge(flag)
       
            
           if(flag == 1)  then
                write(*,*)  "...Converge..."
                exit 
           end if  
           
    end do
    filename ="steady"
    call Output(filename)   
    
    
end subroutine 
