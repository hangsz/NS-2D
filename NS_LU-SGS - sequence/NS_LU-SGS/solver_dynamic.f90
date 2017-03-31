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
    real(8)::oldRsi(5)
    !----------------------------------------------------
    write(*,*)  "solver_dynamic"
    
    call movingInit
    
    dt_r = 2*pi/stepMax/angularFrequency
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
             
            !-------------------------------
            t_total = t_total + dt_r
            call rotation( t_total,U_Rot,AOA )
            !-------------------------------------
            Wn1 = Wn
            Wn = W
            
            do j =1,5
                QW(j,:) = 2.0 /dt_r *vol*Wn(j,:) -1.0/2/dt_r * vol*Wn1(j,:)
            end do
            
            select case( turModelNum ) 
            case( 0 ) 
            case( 1 ) 
                W_SA_n1 = W_SA_n                          ! the n step conservative variables of turbulent control equations       
                W_SA_n = W_SA 
                QW_SA = 2.0 /dt_r *vol*W_SA_n -1.0/2/dt_r * vol*W_SA_n1
            case( 2 )
                W_SST_n1 = W_SST_n                          ! the n step conservative variables of turbulent control equations       
                W_SST_n = W_SST 
                do j=1,2
                    QW_SST(j,:) = 2.0 /dt_r *vol*W_SST_n(j,:) -1.0/2/dt_r * vol*W_SST_n1(j,:)
                end do  
            case default 
            end select  
    
            do count =1,iter_inner 
                write(*,*)  count , iter, i
                write(*,*) "t:",t_total,"dt:",dt_r
                write(*,*) "AoA:",AoA
          
                call meanEdge       !calculate the primative variables on every edge
                call gradient
                !---------------------------------------------------
                ! turbulent residual
                select case( turModelNum ) 
                case( 0 )       
                case( 1 )
                    call getMuT_SA( U_av,muT,muT_av,U_Rot)   ! get turbulent viscosity on the edge
                    call Resi_SA( U,U_av,muL,muL_av,Grad,U_Rot ) 
                    Rsi_SA = Rsi_SA + 3.0/2/dt_r *vol*W_SA - QW_SA
                    call LU_SGS_SA
                case( 2 )
                    call getMUT_SST(U_av,muT,muT_av,U_Rot)
                    call Resi_SST(u_inf,v_inf,U,U_av,muL,muL_av,muT,muT_av,Grad, U_Rot )     !calculate the muT
                    do j=1,2
                        Rsi_SST(j,:) = Rsi_SST(j,:) + 3.0/2/dt_r *vol*W_SST(j,:) - QW_SST(j,:)
                    end do
                    call LU_SGS_SST
                case default 
                end select       
                      
                !-------------------------------------------
                ! renew the turbulent viscosity on the edge
                select case( turModelNum ) 
                case( 1 )
                    call getMuT_SA( U_av,muT,muT_av,U_Rot)
                case( 2 )
                    call getMUT_SST(U_av,muT,muT_av,U_Rot) 
                case default 
            
                end select
                !-------------------------------------------
                ! convective flux
                select case( conFluxNum )
                case( 1 )
                     call conFlux_JST    
                     call artDissipation 
                        
                case( 2 ) 
                     !calculate the convective flux using Roe scheme 
                     Dissi =0.0
                     call conFlux_Roe
                     call timeStep
                end select
                 
                ! viscous flux
        
                call visFlux        !calculate the viscous flux
                !--------------------------
                
                Rsi = Fc - Dissi - Fv
        
                do j =1,5
                    Rsi(j,:) = Rsi(j,:) + 3.0/2/dt_r *vol*W(j,:) - QW(j,:)   
                end do
                 
                call LU_SGS
                !---------------------------------------
          
                if( count==1)  oldRsi = sum( abs(Rsi(1,:)))/ncells
        
                call Converge(flag,oldRsi)
        
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