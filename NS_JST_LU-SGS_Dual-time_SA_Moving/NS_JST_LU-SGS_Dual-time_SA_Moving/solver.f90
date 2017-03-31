subroutine solver          
    !-----------------------------------------------------
        ! purpose: the main function and start funcion 
    !----------------------------------------------------
        
    implicit none
    integer::i,j
    integer::iter                    !iterative variable
    integer::count 
    integer::flag                    !the variable to judge wheathe the mean density converges
    character(len=30)::filename
    real(8)::t_total = 0.0
    real(8)::omg
    real(8)::AoA,angular_rate       ! attact angel, angular velocity
    
    !----------------------------------------------------
    
    write(*,*)  "solver"
     
    call grid                       ! get information about grid              
    call allocateMemory             ! allocate memory for variables
    call allocateMemory_tur
    
    call readFlow                   ! read steady flow properties
    call flowInit                   ! initialize the conservative variables
    
    call reorder                    ! reordering the cell number sequence which is needed for LU-SGS
    
    !--------------------------------------------
    
              ! grid information 
   
    ! store the initial grid information
    vector0 = vector
    rij0 =rij
    rR0 = rR
    rL0 = rL
    tij0 = tij
    
    do i=1,nnodes
        xy0(:,i) = xy(:,i)-rot_cen(:)    ! translate the grid node's coordinates to where the rotational center's coordinates are zero.
    end do 
    
    !-------------------------------------------------
    
    ! angular frequency
     omg =  2.0*kr*Ma_inf*a_inf
     
    !-------------------------------------------------  
     
                         ! dual time scheme
                         
    ! initialize Wn and Wn1
    Wn = W                          ! the n step conservative variables of main control equations
    Wn1 = W                         ! the n+1 step conservative variables of main control equations
    
    WTn=WT                          ! the n step conservative variables of turbulent control equations       
    WTn1=WT                         ! the n+1 step conservative variables of turbulent control equations
    
    !-----------------------------
    do i=1,phase
    
        do iter = 1,itermax
             
            t_total = t_total + dt_r
            AoA = att + att_ampl* sin( omg * t_total )
            angular_rate = omg * ( att_ampl/180.0 *pi ) *cos( omg * t_total)
            
            Wn1 = Wn
            Wn = W
            
            WTn1=WTn
            WTn=WT
                
            do j =1,5
                Q(j,:) = 2.0 /dt_r *vol*Wn(j,:) -1.0/2/dt_r * vol*Wn1(j,:)
            end do
                
            Q_WT = 2.0 /dt_r *vol*WTn -1.0/2/dt_r * vol*WTn1
            
            call Rotation(AoA,angular_rate)
             
            do count =1,iter_inner 
                write(*,*)  count , iter, i
                write(*,*) "t:",t_total,"dt:",dt_r
                write(*,*) "AoA:",AoA
          
                call meanEdge       !calculate the primative variables on every edge
              
                !---------------------------------------------------
                
                ! the viscosity on the edges
                muL =   1.45*( U(5,:)/R/U(1,:) )**(3.0/2)/ ( U(5,:)/R/U(1,:) + 110.0 ) * 1.0E-6
                muL_av = 1.45*U_av(6,:)**(3.0/2)/( U_av(6,:)+110.0 ) *1.0E-6
                
                !turbulence
                call getmuT
                
                !---------------------------------------------------
                
                call conFlux       !calculate the convective flux using Roe scheme
        
                call artDissipation 
        
                dt = CFL*vol/(lamda_c + C_time*lamda_v )        !calculate the time step
        
                call gradient        !calculate the gradient of every cell using the Mean_edge subroutine
           
                call SA_solver       !calculate the muT
                
                !--------------------------
                
                
                ! the residual of turbulent control equations
                RsiT = RsiT + 3.0/2/dt_r *vol*WT - Q_WT
                
                RSiT = 0.0
                call visFlux        !calculate the viscous flux
        
                Rsi = Fc - Dissi - Fv
             
                do j =1,5
                   
                    Rsi(j,:) = Rsi(j,:) + 3.0/2/dt_r *vol*W(j,:) - Q(j,:)
                    
                end do
                
               !-----------------------------------
                
                call LU_SGS
                
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
        call Output(filename) 
        write(*,*) "Output" 
        
    end do
    call outputFreeFlow
    
end subroutine 
