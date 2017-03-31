subroutine  muT_cal
    implicit none
    integer::i
    integer::iter
    integer::flag
    integer::ncl,ncr
    real(8)::Vn
    
    real(8),allocatable::chi_max(:)
   
    allocate( chi_max(ncells) )
    
    chi = nuT/( muL/U(1,:) ) 
    
    muT =  chi**3/( chi**3 + Cv1**3 )* U(1,:)*nuT
     
    do i=1,nedges
        ncl = iedge(3,i)                      !atain the number of the left and the right cell
        ncr = iedge(4,i)
        select case( ncr ) 
        case (-1)                           !-1 represents that the edge is the aerofoil  surface         
            muT_av(i) = 0.0
        case (-2)   !-2 represents that the edge is the farfield  boundaries
            
            Vn = dot_product( U_av(2:3,i) , vector(:,i) ) 
            if ( Vn .LE. 0.0) then
                muT_av(i)=0.1**3/(0.1**3 + Cv1**3) * U_av(1,i)*( 0.1*muL_inf/rou_inf  )  !chi_w = 0.1 = 0.1 *nuL_inf/nuL_inf
            else
                muT_av(i) = muT(ncl)
            end if
            
        case default
        
             muT_av(i)=0.5*( muT(ncl) + muT(ncr) )  
            
        end select
     
    end do
    
    
    
    end subroutine