subroutine  getmuT
    !------------------------------------------------------------
        ! purpose: get the turbulent viscosity coefficient
    !-----------------------------------------------------------
    implicit none
    integer::i
    integer::iter
    integer::flag
    integer::ncl,ncr
    real(8)::ux,uy,Vn
    
    real(8),allocatable::chi_max(:)
   
    allocate( chi_max(ncells) )
    
    chi = nuT/( muL/U(1,:) ) 
    
    muT =  chi**3/( chi**3 + Cv1**3 )* U(1,:)*nuT
     
    do i=1,nedges
        ncl = iedge(3,i)                    
        ncr = iedge(4,i)
        
        select case( ncr ) 
            
        case (-1)    ! -1 represents that the edge is the wall, using the wall boundary condition
            muT_av(i) = 0.0
            
        case (-2)   ! -2 represents that the edge is the farfield  boundary
            
            ux =  U_av(2,i) - 0.5 * ( U_Rot(1,iedge(1,i)) + U_Rot(1,iedge(2,i)) )
            uy =  U_av(3,i) - 0.5 * ( U_Rot(2,iedge(1,i)) + U_Rot(2,iedge(2,i)) )
            Vn = ux* vector(1,i)   + uy * vector(2,i)
            ! Vn = dot_product( U_av(2:3,i) , vector(:,i) ) 
            
            if ( Vn .LE. 0.0) then          ! inflow
                ! using the inflow boundary contion
                muT_av(i)=0.1**3/(0.1**3 + Cv1**3) * U_av(1,i)*( 0.1*muL_inf/rou_inf  )  !chi_w = 0.1 = 0.1 *nuL_inf/nuL_inf
            else                            ! outflow
                ! using the outflow boundary  conditon 
                muT_av(i) = muT(ncl)
            end if
            
        case default   ! edges in the interal domain
        
             muT_av(i)=0.5*( muT(ncl) + muT(ncr) )  
            
        end select
     
    end do
    
    end subroutine