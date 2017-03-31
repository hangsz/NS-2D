subroutine LU_SGS   !Lower-Upper Symmetric Successive Overrelexation scheme
    !object: get W(n+1)  at every time step
    !dual time shceme needs to change the diagonal matrix
    implicit none
    integer::i,j
    
    integer(8)::edge,cell_num,cell_lower,cell_upper
    
    real(8)::sign
    real(8)::Vn,c,rA_v,rA
    real(8)::U_temp(5),W_temp(5),Fc_n1(5),Fc_n(5)
    real(8)::dFc(5)
    real(8)::diag(5)
    real(8),allocatable::dw1(:,:)
    real(8),allocatable::dwn(:,:)
    
    !turbulence
    real(8)::WT_temp(2),Fc_WT_n1(2),Fc_WT_n(2)
    real(8)::dFc_WT(2)
    real(8)::diag_WT(2)
    real(8),allocatable::dWT1(:,:)
    real(8),allocatable::dWTn(:,:)
    
    !turbulent limit
    real(8)::Tu,Tumin,Tumax,kmax,kmin
    real(8)::muT_max,omg_min
    
    !write(*,*)  "LU-SGS_unsteady"
    allocate( dw1(5,ncells) )
    allocate( dwn(5,ncells) )
    
    allocate( dWT1(2,ncells) )
    allocate( dWTn(2,ncells) )
    
    dw1 = 0.0
    dwn = 0.0
    
    dWT1=0.0
    dWTn=0.0
       
    !Forward  sweep
    do i=1,ncells   !it doesn't mean that the sequence is identical to the cell order
        
        cell_num = assigned(i)  
        ! write(*,*)  "Cell assigned:",cell_num , i , ncells
          
        dFc = 0.0  
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
            

            !---------------------------------------------------
            !convective flux
            !Riemann problem's approximate solution, Roe scheme
            
            !calculate  the lower cells's n+1 step convective flux 
            W_temp = W(:,cell_lower) + dw1(:,cell_lower)  
            
            WT_temp = WT(:,cell_lower) + dWT1(:,cell_lower) 
            
            U_temp(1) = W_temp(1)
            U_temp(2) = W_temp(2)/W_temp(1)
            U_temp(3) = W_temp(3)/W_temp(1)
            U_temp(5) = (gamma-1.0)*( W_temp(5)-U_temp(1)*( U_temp(2)**2 + U_temp(3)**2)/2.0 - WT_temp(1))    !E = e + 0.5*V**2 + k
            
          
            Vn = dot_product(U_temp(2:3),vector(:,edge) )
            
            Fc_n1(1) = W_temp(1) * Vn
            Fc_n1(2) = W_temp(2) * Vn + vector(1,edge) * U_temp(5)
            Fc_n1(3) = W_temp(3) * Vn + vector(2,edge) * U_temp(5)
            Fc_n1(5) = ( W_temp(5) + U_temp(5) ) * Vn
            
            Fc_WT_n1 = WT_temp * Vn
            
            !-----------------------------
            !calculate  the lower cells's n step convective flux
         
            Vn = dot_product( U(2:3,cell_lower),vector(:,edge) )
            c = sqrt( gamma*U(5,cell_lower) / U(1,cell_lower) )
            
            rA_v = ds(edge) /sqrt( sum( (cell_center(:,cell_num)-cell_center(:,cell_lower) )**2) )* max( 4.0/3/U(1,cell_lower) ,gamma/U(1,cell_lower) )*( muL(cell_lower)/PrL + muT(cell_lower)/PrT ) 
            rA = omg * ( abs(Vn) + c*ds(edge) ) + rA_v
        
            Fc_n(1) = W(1,cell_lower) * Vn
            Fc_n(2) = W(2,cell_lower) * Vn + vector(1,edge) * U(5,cell_lower)
            Fc_n(3) = W(3,cell_lower) * Vn + vector(2,edge) * U(5,cell_lower)
            Fc_n(5) = ( W(5,cell_lower) + U(5,cell_lower) )* Vn
            
            Fc_WT_n = WT(:,cell_lower) *Vn
            
            dFc = dFc + sign*(Fc_n1 - Fc_n ) -  rA*dw1(:,cell_lower)  
            
            dFc_WT = dFc_WT + sign*(Fc_WT_n1 - Fc_WT_n ) -  rA*dWT1(:,cell_lower)     
            
        end do 
        
        diag = vol(cell_num)/dt(cell_num)  + omg/2.0*lamda_c(cell_num) + lamda_v(cell_num)
        
        diag_WT = vol(cell_num)/dt(cell_num)  + omg/2.0*lamda_c(cell_num) + lamda_v(cell_num)  - vol(cell_num)*dQTdWT(:,cell_num) 
        
        dw1(:,cell_num) = 1.0/diag * ( - Rsi(:,cell_num) - 0.5*dFc )
        
        dWT1(:,cell_num) = 1.0/diag_WT * ( - RsiT(:,cell_num) - 0.5*dFc_WT )

    end do    
    
    !Backward sweep
    do i = ncells,1,-1
        
        cell_num = assigned(i)    
        
        dFc = 0.0
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
                         
            !-----------------------------------------------
            !calculate the upper cells's n+1 convective flux 
            W_temp = W(:,cell_upper) + dwn(:,cell_upper)  
            WT_temp = WT(:,cell_upper) + dWTn(:,cell_upper)  
            
            U_temp(1) = W_temp(1)
            U_temp(2) = W_temp(2)/W_temp(1)
            U_temp(3) = W_temp(3)/W_temp(1)
            U_temp(5) = (gamma-1.0)*( W_temp(5)-U_temp(1)*( U_temp(2)**2 + U_temp(3)**2)/2.0 - WT_temp(1) )                            !E = e + 0.5*V**2 + k
            
            Vn = dot_product( U_temp(2:3),vector(:,edge) )
            
            Fc_n1(1) = W_temp(1) * Vn
            Fc_n1(2) = W_temp(2) * Vn + vector(1,edge) *U_temp(5)
            Fc_n1(3) = W_temp(3) * Vn + vector(2,edge) *U_temp(5)
            Fc_n1(5) = ( W_temp(5) + U_temp(5) ) * Vn
         
            Fc_WT_n1 = WT_temp * Vn
            
            !------------------------------
            !calculate  the upper cells's n step convective flux 
            
            Vn = dot_product( U(2:3,cell_upper),vector(:,edge) )
            c = sqrt( gamma*U(5,cell_upper) / U(1,cell_upper) )
            rA_v = ds(edge) / sqrt( sum( (cell_center(:,cell_num)-cell_center(:,cell_upper) )**2 ) )*max( 4.0/3/U(1,cell_upper) ,gamma/U(1,cell_upper) )*( muL(cell_upper)/PrL + muT(cell_upper)/PrT ) 
            rA =omg*(  abs(Vn) + c*ds(edge) ) + rA_v
            
            
            Fc_n(1) = W(1,cell_upper) * Vn
            Fc_n(2) = W(2,cell_upper) * Vn + vector(1,edge) *U(5,cell_upper)
            Fc_n(3) = W(3,cell_upper) * Vn + vector(2,edge) *U(5,cell_upper)
            Fc_n(5) = ( W(5,cell_upper) + U(5,cell_upper) )* Vn
            
            Fc_WT_n = WT(:,cell_upper) *Vn
            
            dFc = dFc + sign*( Fc_n1 - Fc_n ) - rA*dwn(:,cell_upper)        
            dFc_WT = dFc_WT + sign*( Fc_WT_n1 - Fc_WT_n ) - rA*dWTn(:,cell_upper)          
            
        end do
        
        diag = vol(cell_num)/dt(cell_num)  + omg/2.0*lamda_c(cell_num) + lamda_v(cell_num)
        
        diag_WT = vol(cell_num)/dt(cell_num)  + omg/2.0*lamda_c(cell_num) + lamda_v(cell_num)  - vol(cell_num)*dQTdWT(:,cell_num) 
        
        dwn(:,cell_num) = dw1(:,cell_num) -  0.5*dFc/diag
        
        dWTn(:,cell_num) = dWT1(:,cell_num) -  0.5*dFc_WT/diag_WT
        
    end do
    
    !calculate the original variables
    W = W + dwn 
    WT = WT + dWTn  
    
    U(1,:) = W(1,:)
    U(2,:) = W(2,:)/W(1,:)
    U(3,:) = W(3,:)/W(1,:)
    !U(5,:) = (gamma-1.0)*( W(5,:)-U(1,:)*( U(2,:)**2 + U(3,:)**2)/2.0 -WT(1,:) )  
    
    !do i=1,ncells
      
        !Tu = 0.16* ( U(1,i)*sqrt( U(2,i)**2 + U(3,i)**2 ) /muT(i)  )**(-1.0/8)
        !Tumin = min( 1.0E-8 , 0.001*Tu )
        !Tumax = max( 0.8 , 100.0*Tu)
        !kmin = 1.5 * Tumin**2*( U(2,i)**2 + U(3,i)**2)
        !kmax = 1.5 * Tumax**2*( U(2,i)**2 + U(3,i)**2 + 0.01*Ma_inf*a_inf )
         
        
        !UT(1,i) = WT(1,i) / U(1,i)
        
        !if ( UT(1,i) .GT. kmax )  then
        !    UT(1,i) = kmax
        !end if
        !if ( UT(1,i) .LT. kmin ) then
        !    UT(1,i) = kmin
        !end if
        
        !muT_max = 5000.0 * muL(i)
        !omg_min = max( sqrt( UT(1,i)/d(i)**2 ),U(1,i)*UT(1,i)/muT_max )
        !!
        !UT(2,i) = WT(2,i) / U(1,i)
        !if ( UT(2,i) .LT. omg_min ) then
        !    UT(2,i) = omg_min
        !end if
        !
        !
        !if( UT(1,i)<0.0 .OR. UT(2,i)<0.0 ) then
        !    write(*,*) UT(:,i)
        !    stop
        !end if
              
    !end do
    
    UT(1,:) = WT(1,:) / U(1,:)
    UT(2,:) = WT(2,:) / U(1,:)
    
    where( UT(1,:)<= 0.0 )
        UT(1,:) = minval(abs(UT(1,:)))
    end where
        
    where( UT(2,:)<= 0.0 )
        UT(2,:) = minval(abs(UT(2,:)))
    end where
    
    do i=1, ncells  
        muT(i) = a1*U(1,i)*UT(1,i)/max( a1*UT(2,i),f2(i)*abs( Grad(2,2,i) - Grad(1,3,i)) )  
    end do
    
    U(5,:) = (gamma-1.0)*( W(5,:)-U(1,:)*( U(2,:)**2 + U(3,:)**2)/2.0 -U(1,:)*UT(1,:) )  
    
end subroutine
        
        
        
        
        