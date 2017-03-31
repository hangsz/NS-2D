subroutine solver_dynamic          
    !-----------------------------------------------------
        ! purpose: the main function and start funcion 
    !----------------------------------------------------
        
    implicit none
    integer::i,j
    integer::iter                    !iterative variable
    integer::count 
    integer::flag                    !the variable to judge wheathe the mean density converges
    character(len=50)::filename
    real(8)::t_total = 0.0
    real(8)::AoA,angular_rate       ! attact angel, angular velocity
    
    !----------------------------------------------------
    write(*,*)  "solver_dynamic"
    
    call movingInit
    
    
    dt_r = pi/2/stepMax/angularFrequency
    

    ! dual time                     
    ! initialize Wn and Wn1
    Wn = W                                          ! the n step conservative variables of main control equations
    Wn1 = W                                         ! the n+1 step conservative variables of main control equations
    
    !------------------------------
    select case( turModelNum )      
    case( 0 )        
        case( 1 ) 
            W_SA_n = W_SA                            ! the n step conservative variables of turbulent control equations       
            W_SA_n1 = W_SA   
        case( 2 )
            W_SST_n = W_SST                          ! the n step conservative variables of turbulent control equations       
            W_SST_n1 = W_SST  
        case default 
    end select        
        
    !-----------------------------
    do i=1,phase
    
        do iter = 1,steps
             
            t_total = t_total + dt_r
            AoA = att + att_ampl* sin( angularFrequency * t_total )
            angular_rate = angularFrequency * ( att_ampl/180.0 *pi ) *cos( angularFrequency * t_total)
           
            !-------------------------------------
            Wn1 = Wn
            Wn = W
            select case( turModelNum ) 
                case( 0 ) 
                case( 1 ) 
                    W_SA_n1 = W_SA_n                          ! the n step conservative variables of turbulent control equations       
                    W_SA_n = W_SA   
                case( 2 )
                    W_SST_n1 = W_SST_n                          ! the n step conservative variables of turbulent control equations       
                    W_SST_n = W_SST  
                 case default 
            end select  
            
            !------------------------------------
            do j =1,5
                QW(j,:) = 2.0 /dt_r *vol*Wn(j,:) -1.0/2/dt_r * vol*Wn1(j,:)
            end do
            
            
            select case( turModelNum ) 
                case( 0 ) 
                case( 1 ) 
                   QW_SA = 2.0 /dt_r *vol*W_SA_n -1.0/2/dt_r * vol*W_SA_n1
                case( 2 )
                    do j=1,2
                       QW_SST(j,:) = 2.0 /dt_r *vol*W_SST_n(j,:) -1.0/2/dt_r * vol*W_SST_n1(j,:)
                    end do  
                 case default 
            end select  
                
           !-----------------------------------
            call rotation( t_total,U_Rot )
             
            do count =1,iter_inner 
                write(*,*)  count , iter, i
                write(*,*) "t:",t_total,"dt:",dt_r
                write(*,*) "AoA:",AoA
          
                call meanEdge       !calculate the primative variables on every edge
              
                !---------------------------------------------------
                ! the viscosity on the edges
                muL =   1.45*( U(5,:)/R/U(1,:) )**(3.0/2)/ ( U(5,:)/R/U(1,:) + 110.0 ) * 1.0E-6
                muL_av = 1.45*U_av(6,:)**(3.0/2)/( U_av(6,:)+110.0 ) *1.0E-6
                
                !------------------------------------------
                ! calculate the turbulent properies on the edge
                select case( turModelNum ) 
                case( 1 )
                    
                    call getMuT_SA( U,U_av,muL,muT,muT_av,U_Rot)
                    call meanEdge_SA( U_av,U_Rot)  
                case( 2 )
                    
                    call getMUT_SST(U,U_av,muL,muT,muT_av,U_Rot)
                    call meanEdge_SST( u_inf,v_inf, U_av,muL_av ,U_Rot)   ! its needed to calculate convective flux 
                case default 
            
                end select
                !-------------------------------------------
                ! convective flux
                select case( conFluxNum )
                case( 1 )
                     call conFlux_JST    
         
                     call artDissipation 
                     dt = CFL*vol/(lamda_c + C_time*lamda_v )     
                case( 2 ) 
                     !calculate the convective flux using Roe scheme 
                     call conFlux_Roe
                     call timeStep
                end select
                 
                !--------------------------
                ! viscous flux
                call gradient        !calculate the gradient of every cell using the Mean_edge subroutine 
        
                call visFlux        !calculate the viscous flux
                !--------------------------
                
                Rsi = Fc - Dissi - Fv
        
                do j =1,5
                    Rsi(j,:) = Rsi(j,:) + 3.0/2/dt_r *vol*W(j,:) - QW(j,:)   
                end do
                 
                !---------------------------------------
                ! turbulent residual
                select case( turModelNum ) 
                case( 0 ) 
                         call LU_SGS
                    case( 1 )
                         call Resi_SA( U,U_av, alf,mul,mul_av ,Grad , U_Rot ) 
                         Rsi_SA = Rsi_SA + 3.0/2/dt_r *vol*W_SA - QW_SA
                         call LU_SGS_SA
                    case( 2 )
                         call Resi_SST( U,U_av, alf,muL,muL_av,muT,muT_av,Grad, U_Rot )     !calculate the muT
                         do j=1,2
                             Rsi_SST(j,:) = Rsi_SST(j,:) + 3.0/2/dt_r *vol*W_SST(j,:) - QW_SST(j,:)
                         end do
                         call LU_SGS_SST
                    case default 
                end select       
          
                if( count==1)  oldRsi = sum( abs(Rsi(1,:)))/ncells
        
                call Converge(flag)
        
                write(*,*)
                if(flag == 1)  then
                    write(*,*)  "...Converge..."
                    exit 
                end if 
                
            end do
        end do
        
        write(filename,"(I2)") i
            
        filename = "flow-info-" //trim(filename)
        call output(filename)    
    end do
    
end subroutine 