subroutine LU_SGS_SA!( W, U, assigned, lower,upper)   !Lower-Upper Symmetric Successive Overrelexation scheme
    !object: get W(n+1)  at every time step
    !dual time shceme needs to change the diagonal matrix
    implicit none
    
    !integer,intent(in)::W(:,:),U(:,:), assigned(:),lower(:,:),upper(:,:)
    integer::i,j
    
    integer(8)::edge,cell_num,cell_lower,cell_upper
    
    real(8)::sign
    real(8)::ut,vt
    real(8)::ux,uy,Vn,c,rA_v,rA
    
    !turbulence
    real(8)::WT_temp,Fc_WT_n1,Fc_WT_n
    real(8)::dFc_WT
    real(8)::diag_WT
    real(8),allocatable::dWT1(:)
    real(8),allocatable::dWTn(:)
    
    !write(*,*)  "LU-SGS_unsteady"
    
    allocate( dWT1(ncells) )
    allocate( dWTn(ncells) )
    
    dWT1=0.0
    dWTn=0.0
       
    !Forward  sweep
    do i=1,ncells   !it doesn't mean that the sequence is identical to the cell order
        
        cell_num = assigned(i)  
        ! write(*,*)  "Cell assigned:",cell_num , i , ncells
          
        dFc_WT = 0.0
         
        do j=1,4 
            edge = lower( j,cell_num )   !each cell has three lower edges totally(in fact,it can only have two lower or upper edges at most)   
            if( edge .EQ. 0)  then
            !    write(*,*) "edge:",edge,j
                exit
            end if 
            !determine the cell number lower than the current cell
            !because only the lower edge is stored , it has to be tested which cell is the lower cell
            ! write(*,*)  "ncl:",iedge(3,edge), "ncr:",iedge(4,edge),"cell_num:",cell_num
            if( iedge(3,edge) .NE. cell_num ) then
                cell_lower = iedge(3,edge)
                sign =  -1.0
            !    write(*,*)  "ncl"
            else
                cell_lower = iedge(4,edge)
                !write(*,*)  "ncr"
                sign = 1.0
            end if
            !-------------------------------------
            !convective flux , calculate  the lower cells's n+1 step convective flux 
           
            WT_temp = W_SA(cell_lower) + dWT1(cell_lower)  
            
            if( icell(3,cell_lower)/=icell(4,cell_lower) ) then
                ut = 1.0/4 * sum( U_Rot(1,icell(:,cell_lower)) )
                vt = 1.0/4 * sum( U_Rot(2,icell(:,cell_lower)) )
            else
                ut = 1.0/3 * sum( U_Rot(1,icell(1:3,cell_lower)) )
                vt = 1.0/3 * sum( U_Rot(2,icell(1:3,cell_lower)) )
            end if
            
            ux =  U(2,cell_lower) - ut
            uy =  U(3,cell_lower) - vt
            Vn = ux * vector(1,edge)   + uy * vector(2,edge) 
    
            Fc_WT_n1 = WT_temp * Vn
            !------------------
            !calculate  the lower cells'S_SA n step convective flux 
           
            rA = omg*alf_SA(edge) 
            Fc_WT_n = W_SA(cell_lower) *Vn
            
            dFc_WT = dFc_WT + sign*(Fc_WT_n1 - Fc_WT_n ) -  rA*dWT1(cell_lower)     
            
        end do 
        
        !-------------------------------------------------------------------
    
        if( moveNum == 0 ) then
            diag_WT = vol(cell_num)/dt(cell_num)  + omg/2.0*lamda_SA(cell_num) - vol(cell_num)*dQdW_SA(cell_num) 
        else
            diag_WT = vol(cell_num)*( 1.0/dt(cell_num) + 3.0/2/dt_r ) + omg/2.0*lamda_SA(cell_num)  - vol(cell_num)*dQdW_SA(cell_num) 
        end if
        
        dWT1(cell_num) = 1.0/diag_WT * ( - Rsi_SA(cell_num) - 0.5*dFc_WT )

    end do    
    
    !Backward sweep
    do i = ncells,1,-1
        
        cell_num = assigned(i)    
        
        dFc_WT=0.0
        
        do j=1,4
            
            edge = upper( j,cell_num )
            
            if( edge .EQ. 0 )  exit
                
            if( iedge(3,edge) .NE. cell_num ) then
                cell_upper = iedge(3,edge)
                sign = -1.0
            else
                cell_upper = iedge(4,edge)
                sign = 1.0
            end if
                         
            !-------------------------------------------------
            !calculate the upper cells's n+1 convective flux 
            WT_temp = W_SA(cell_upper) + dWTn(cell_upper)  
            
            if( icell(3,cell_upper)/=icell(4,cell_upper) ) then
                ut = 1.0/4 * sum( U_Rot(1,icell(:,cell_upper)) )
                vt = 1.0/4 * sum( U_Rot(2,icell(:,cell_upper)) )
            else
                ut = 1.0/3 * sum( U_Rot(1,icell(1:3,cell_upper)) )
                vt = 1.0/3 * sum( U_Rot(2,icell(1:3,cell_upper)) )
            end if
            
            ux =  U(2,cell_upper) - ut
            uy =  U(3,cell_upper) - vt
            Vn = ux * vector(1,edge)   + uy * vector(2,edge) 
         
            Fc_WT_n1 = WT_temp * Vn
            
            !calculate  the upper cells's n step convective flux 
            
            rA= omg*alf_SA(edge) 
            Fc_WT_n = W_SA(cell_upper) *Vn
            
            dFc_WT = dFc_WT + sign*( Fc_WT_n1 - Fc_WT_n ) - rA*dWTn(cell_upper)          
            
        end do
        
        if( moveNum == 0 ) then
            diag_WT = vol(cell_num)/dt(cell_num)  + omg/2.0*lamda_SA(cell_num)  - vol(cell_num)*dQdW_SA(cell_num) 
        else
            diag_WT = vol(cell_num)*( 1.0/dt(cell_num) + 3.0/2/dt_r ) + omg/2.0*lamda_SA(cell_num)  - vol(cell_num)*dQdW_SA(cell_num) 
        end if
        dWTn(cell_num) = dWT1(cell_num) -  0.5*dFc_WT/diag_WT
        
    end do
    
    W_SA = W_SA + dWTn 
        
    U_SA = W_SA/U(1,:)
    where( U_SA <= 0.0 ) 
        U_SA = minval( abs(U_SA))
    end where
    
    where( U_SA >( 5000.0*muL_inf/rou_inf ) )
        U_SA = 5000.0*muL_inf/rou_inf
    end where
    
    chi_SA = U_SA/( muL/U(1,:) ) 
    muT =  chi_SA**3/( chi_SA**3 + Cv1_SA**3 )* U(1,:)*U_SA
    
end subroutine
        
        
        
        
        