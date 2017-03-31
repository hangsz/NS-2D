subroutine LU_SGS   !Lower-Upper Symmetric Successive Overrelexation scheme
    !object: get W(n+1)  at every time step
    !dual time shceme needs to change the diagonal matrix
    implicit none
    integer::i,j
    
    integer(8)::edge,cell_num,cell_lower,cell_upper
    
    real(8)::sign
    real(8)::ut,vt,Vnt
    real(8)::ux,uy,Vn,c,rA,rA_V
    real(8)::U_temp(5),W_temp(5),Fc_n1(5),Fc_n(5)
    real(8)::dFc(5)
    
    real(8)::diag(5)
    real(8)::drou_old = 0.0 ,drou_new
    
    
    real(8),allocatable::dw1(:,:)
    real(8),allocatable::dwn(:,:)
    
    !write(*,*)  "LU-SGS_unsteady"
    allocate( dw1(5,ncells) )
    allocate( dwn(5,ncells) )
    
    dw1 = 0.0
    dwn = 0.0
       
    !Forward  sweep
    do i=1,ncells   !it doesn't mean that the sequence is identical to the cell order
        
        cell_num = assigned(i)  
        ! write(*,*)  "Cell assigned:",cell_num , i , ncells
          
        dFc = 0.0  
         
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
            
            !----------------------------------------------------------------
            !convective flux
            !---------------------------------
            !calculate  the lower cells's n+1 step convective flux 
            W_temp = W(:,cell_lower) + dw1(:,cell_lower)  
         
            U_temp(1) = W_temp(1)
            U_temp(2) = W_temp(2)/W_temp(1)
            U_temp(3) = W_temp(3)/W_temp(1)
            U_temp(5) = (gamma-1.0)*( W_temp(5)-U_temp(1)*( U_temp(2)**2 + U_temp(3)**2)/2.0 ) 
           
            if( icell(3,cell_lower)/=icell(4,cell_lower) ) then
                ut = 1.0/4 * sum( U_Rot(1,icell(:,cell_lower)) )
                vt = 1.0/4 * sum( U_Rot(2,icell(:,cell_lower)) )
            else
                ut = 1.0/3 * sum( U_Rot(1,icell(1:3,cell_lower)) )
                vt = 1.0/3 * sum( U_Rot(2,icell(1:3,cell_lower)) )
            end if
            
            Vnt = ut*vector(1,edge) + vt*vector(2,edge)
            
            ux =  U_temp(2) - ut
            uy =  U_temp(3) - vt
            !ux =  U_temp(2) - 1.0/4 * sum( U_Rot(1,icell(:,cell_lower)) )
            !uy =  U_temp(3) - 1.0/4 * sum( U_Rot(2,icell(:,cell_lower)) )
            Vn = ux * vector(1,edge)   + uy * vector(2,edge) 
            !Vn = dot_product(U_temp(2:3),vector(:,edge) )
            
            Fc_n1(1) = W_temp(1) * Vn
            Fc_n1(2) = W_temp(2) * Vn + vector(1,edge) * U_temp(5)
            Fc_n1(3) = W_temp(3) * Vn + vector(2,edge) * U_temp(5)
            Fc_n1(5) = ( W_temp(5) + U_temp(5) ) * Vn
            Fc_n1(5) = Fc_n1(5) + Vnt*U_temp(5)
            !--------------------
            !calculate  the lower cells's n step convective flux 
            ux =  U(2,cell_lower) - ut
            uy =  U(3,cell_lower) - vt
            !ux =  U(2,cell_lower) - 1.0/4 * sum( U_Rot(1,icell(:,cell_lower)) )
            !uy =  U(3,cell_lower) - 1.0/4 * sum( U_Rot(2,icell(:,cell_lower)) )
            Vn = ux * vector(1,edge)   + uy * vector(2,edge)
            !Vn = dot_product( U(2:3,cell_lower),vector(:,edge) )
            c = sqrt( gamma*U(5,cell_lower) / U(1,cell_lower) )
         
            rA_v = ds(edge) /sqrt( ( cell_center(1,cell_num)-cell_center(1,cell_lower) )**2 + &
                 & ( cell_center(2,cell_num)-cell_center(2,cell_lower) )**2 )* max( 4.0/3/U(1,cell_lower) ,gamma/U(1,cell_lower) )*( muL(cell_lower)/PrL + muT(cell_lower)/PrT ) 
            rA = omg * ( abs(Vn) + c*ds(edge) )  +rA_v
           
            Fc_n(1) = W(1,cell_lower) * Vn
            Fc_n(2) = W(2,cell_lower) * Vn + vector(1,edge) * U(5,cell_lower)
            Fc_n(3) = W(3,cell_lower) * Vn + vector(2,edge) * U(5,cell_lower)
            Fc_n(5) = ( W(5,cell_lower) + U(5,cell_lower) )* Vn
            Fc_n(5) = Fc_n(5) + Vnt * U(5,cell_lower)
            
            dFc = dFc + sign*(Fc_n1 - Fc_n ) -  rA*dw1(:,cell_lower)    
        end do 
        
        if( moveNum==0 )  then
            diag = vol(cell_num)/dt(cell_num)+ omg/2.0*lamda_c(cell_num) + lamda_v(cell_num)
        else
            diag = vol(cell_num)*(1.0/dt(cell_num) + 3.0/2/dt_r ) + omg/2.0*lamda_c(cell_num) + lamda_v(cell_num)
        end if
        !diag = vol(cell_num)*( 1.0/dt(cell_num) + 3.0/2/dt_r ) + omg/2.0*lamda_c(cell_num) + lamda_v(cell_num) 
        
        dw1(:,cell_num) = 1.0/diag * ( - Rsi(:,cell_num) - 0.5*dFc )

    end do    
    
    !Backward sweep
    do i = ncells,1,-1
        
        cell_num = assigned(i)    
        
        dFc = 0.0

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
                         
            !---------------------------------------------------------------------
            !calculate the upper cells's n+1 convective flux 
            W_temp = W(:,cell_upper) + dwn(:,cell_upper)  
         
            U_temp(1) = W_temp(1)
            U_temp(2) = W_temp(2)/W_temp(1)
            U_temp(3) = W_temp(3)/W_temp(1)
            U_temp(5) = (gamma-1.0)*( W_temp(5)-U_temp(1)*( U_temp(2)**2 + U_temp(3)**2)/2.0 )  
           
            if( icell(3,cell_upper)/=icell(4,cell_upper) ) then
                ut = 1.0/4 * sum( U_Rot(1,icell(:,cell_upper)))
                vt = 1.0/4 * sum( U_Rot(2,icell(:,cell_upper)) )
            else
                ut = 1.0/3 * sum( U_Rot(1,icell(1:3,cell_upper)) )
                vt = 1.0/3 * sum( U_Rot(2,icell(1:3,cell_upper)) )
            end if
            
            Vnt = ut * vector(1,edge)   + vt * vector(2,edge)
        
            ux =  U_temp(2) - ut
            uy =  U_temp(3) - vt
            !ux =  U_temp(2) - 1.0/4 * sum( U_Rot(1,icell(:,cell_upper)) )
            !uy =  U_temp(3) - 1.0/4 * sum( U_Rot(2,icell(:,cell_upper)) )    
            Vn = ux * vector(1,edge)   + uy * vector(2,edge)
            !Vn = dot_product( U_temp(2:3),vector(:,edge) )
            
            Fc_n1(1) = W_temp(1) * Vn
            Fc_n1(2) = W_temp(2) * Vn + vector(1,edge) *U_temp(5)
            Fc_n1(3) = W_temp(3) * Vn + vector(2,edge) *U_temp(5)
            Fc_n1(5) = ( W_temp(5) + U_temp(5) ) * Vn
            Fc_n1(5) = Fc_n1(5) + Vnt*U_temp(5)
            !------------------------------------------
            
            !calculate  the upper cells's n step convective flux 
            ux =  U(2,cell_upper) - ut
            uy =  U(3,cell_upper) - vt
            !ux =  U(2,cell_upper) - 1.0/4 * sum( U_Rot(1,icell(:,cell_upper)) )
            !uy =  U(3,cell_upper) - 1.0/4 * sum( U_Rot(2,icell(:,cell_upper)) )   
            Vn = ux * vector(1,edge)   + uy * vector(2,edge)
            !Vn = dot_product( U(2:3,cell_upper),vector(:,edge) )
            c = sqrt( gamma*U(5,cell_upper) / U(1,cell_upper) )
            rA_v = ds(edge) / sqrt( ( cell_center(1,cell_num)-cell_center(1,cell_upper) )**2  +  &
                 & ( cell_center(2,cell_num)-cell_center(2,cell_upper) )**2 )*max( 4.0/3/U(1,cell_upper) ,gamma/U(1,cell_upper) )*( muL(cell_upper)/PrL + muT(cell_upper)/PrT ) 
            rA =omg*(  abs(Vn) + c*ds(edge) ) +rA_v
            
            Fc_n(1) = W(1,cell_upper) * Vn
            Fc_n(2) = W(2,cell_upper) * Vn + vector(1,edge) *U(5,cell_upper)
            Fc_n(3) = W(3,cell_upper) * Vn + vector(2,edge) *U(5,cell_upper)
            Fc_n(5) = ( W(5,cell_upper) + U(5,cell_upper) )* Vn
            Fc_n(5) = Fc_n(5) +  Vnt*U(5,cell_upper) 
            
            dFc = dFc + sign*( Fc_n1 - Fc_n ) - rA*dwn(:,cell_upper)           
            
        end do
        
        if( moveNum==0 )  then
            diag = vol(cell_num)/dt(cell_num)+ omg/2.0*lamda_c(cell_num) + lamda_v(cell_num)
        else
            diag = vol(cell_num)*(1.0/dt(cell_num) + 3.0/2/dt_r ) + omg/2.0*lamda_c(cell_num) + lamda_v(cell_num)
        end if
        
        dwn(:,cell_num) = dw1(:,cell_num) -  0.5*dFc/diag
        
     
    end do
          
    W = W + dwn
    
    !calculate the original variables
    
    U(1,:) = W(1,:)
    U(2,:) = W(2,:)/W(1,:)
    U(3,:) = W(3,:)/W(1,:)
    U(5,:) = (gamma-1.0)*( W(5,:)-U(1,:)*( U(2,:)**2 + U(3,:)**2)/2.0 )  
        
      
end subroutine
        
        
        