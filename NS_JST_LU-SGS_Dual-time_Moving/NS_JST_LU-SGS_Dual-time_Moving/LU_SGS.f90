subroutine LU_SGS   !Lower-Upper Symmetric Successive Overrelexation scheme
    !object: get W(n+1)  at every time step
    !dual time shceme needs to change the diagonal matrix
    implicit none
    integer::i,j
    
    integer(8)::edge,cell_num,cell_lower,cell_upper
    
    real(8)::sign
    real(8)::Vn,c,rA
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
            

            !*****************************************************
            !convective flux
            !Riemann problem's approximate solution, Roe scheme
            
            !********************
            !calculate  the lower cells's n+1 step convective flux 
            W_temp = W(:,cell_lower) + dw1(:,cell_lower)  
         
            U_temp(1) = W_temp(1)
            U_temp(2) = W_temp(2)/W_temp(1)
            U_temp(3) = W_temp(3)/W_temp(1)
            U_temp(5) = (gamma-1.0)*( W_temp(5)-U_temp(1)*( U_temp(2)**2 + U_temp(3)**2)/2.0 )  
            
            Vn = dot_product(U_temp(2:3),vector(:,edge) )
            
            Fc_n1(1) = W_temp(1) * Vn
            Fc_n1(2) = W_temp(2) * Vn + vector(1,edge) * U_temp(5)
            Fc_n1(3) = W_temp(3) * Vn + vector(2,edge) * U_temp(5)
            Fc_n1(5) = ( W_temp(5) + U_temp(5) ) * Vn
            
            !***********
            !calculate  the lower cells's n step convective flux 
            Vn = dot_product( U(2:3,cell_lower),vector(:,edge) )
            c = sqrt( gamma*U(5,cell_lower) / U(1,cell_lower) )
            rA = omg_o * ( abs(Vn) + c*ds(edge) )
           
            
            Fc_n(1) = W(1,cell_lower) * Vn
            Fc_n(2) = W(2,cell_lower) * Vn + vector(1,edge) * U(5,cell_lower)
            Fc_n(3) = W(3,cell_lower) * Vn + vector(2,edge) * U(5,cell_lower)
            Fc_n(5) = ( W(5,cell_lower) + U(5,cell_lower) )* Vn
            
            dFc = dFc + sign*(Fc_n1 - Fc_n ) -  rA*dw1(:,cell_lower)   
               
            
        end do 
        
        
        diag = vol(cell_num)*( 1.0/dt(cell_num) + 3.0/2/dt_r ) + omg_o/2.0*lamda_c(cell_num) + lamda_v(cell_num) 
        
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
                         
            !*****************************************************
            !calculate the upper cells's n+1 convective flux 
            W_temp = W(:,cell_upper) + dwn(:,cell_upper)  
         
            U_temp(1) = W_temp(1)
            U_temp(2) = W_temp(2)/W_temp(1)
            U_temp(3) = W_temp(3)/W_temp(1)
            U_temp(5) = (gamma-1.0)*( W_temp(5)-U_temp(1)*( U_temp(2)**2 + U_temp(3)**2)/2.0 )  
                    
            Vn = dot_product( U_temp(2:3),vector(:,edge) )
            
            Fc_n1(1) = W_temp(1) * Vn
            Fc_n1(2) = W_temp(2) * Vn + vector(1,edge) *U_temp(5)
            Fc_n1(3) = W_temp(3) * Vn + vector(2,edge) *U_temp(5)
            Fc_n1(5) = ( W_temp(5) + U_temp(5) ) * Vn
         
            !***********
            !calculate  the upper cells's n step convective flux 
            Vn = dot_product( U(2:3,cell_upper),vector(:,edge) )
            c = sqrt( gamma*U(5,cell_upper) / U(1,cell_upper) )
            rA =omg_o*(  abs(Vn) + c*ds(edge) )
            
            
            Fc_n(1) = W(1,cell_upper) * Vn
            Fc_n(2) = W(2,cell_upper) * Vn + vector(1,edge) *U(5,cell_upper)
            Fc_n(3) = W(3,cell_upper) * Vn + vector(2,edge) *U(5,cell_upper)
            Fc_n(5) = ( W(5,cell_upper) + U(5,cell_upper) )* Vn
            
            
            dFc = dFc + sign*( Fc_n1 - Fc_n ) - rA*dwn(:,cell_upper)           
            
        end do
        
        diag = vol(cell_num)*(1.0/dt(cell_num) + 3.0/2/dt_r ) + omg_o/2.0*lamda_c(cell_num) + lamda_v(cell_num)
        dwn(:,cell_num) = dw1(:,cell_num) -  0.5*dFc/diag
        
     
    end do
    
           
    W = W + dwn
    
    !calculate the original variables
    
    U(1,:) = W(1,:)
    U(2,:) = W(2,:)/W(1,:)
    U(3,:) = W(3,:)/W(1,:)
    U(5,:) = (gamma-1.0)*( W(5,:)-U(1,:)*( U(2,:)**2 + U(3,:)**2)/2.0 )  
        
      
end subroutine
        
        
        