subroutine Mean_edge_Tur  
    !object: get the flow properties on every edge, it's needed to calculate the gradient
    !        
    implicit none
    integer::i
    integer::ncl,ncr
    real(8)::Vn
    
    nuT_av = 0.0
    
    do i=1,nedges
        ncl = iedge(3,i)                      ! the number of the left and the right cell
        ncr = iedge(4,i)
        select case( ncr ) 
        case (-1)                             !-1 represents that the edge is the aerofoil  surface         
            nuT_av(i) = 0.0
        case (-2)  !-2 represents that the edge is the farfield  boundaries
            
           
            Vn = dot_product( U_av(2:3,i) , vector(:,i) ) 
            if ( Vn .LE. 0.0) then                     
                nuT_av(i)= 0.1 * muL_inf / rou_inf
            else
                nuT_av(i) = nuT(ncl)
            end if
            
        case default
        
            nuT_av(i)=0.5*( nuT(ncl) + nuT(ncr) )  
            
        end select
          
    end do
    
end subroutine