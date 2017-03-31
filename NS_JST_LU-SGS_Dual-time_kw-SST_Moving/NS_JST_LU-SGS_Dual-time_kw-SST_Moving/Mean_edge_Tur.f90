subroutine Mean_edge_Tur  
    !object: get the flow properties on every edge, it's needed to calculate the gradient
    !        
    implicit none
    integer::i
    integer::ncl,ncr
    real(8)::ux,uy,Vn
    
    UT_av = 0.0
    
    do i=1,nedges
        ncl = iedge(3,i)                      ! the number of the left and the right cell
        ncr = iedge(4,i)
        select case( ncr ) 
        case (-1)     !-1 represents that the edge is the aerofoil  surface         
            UT_av(1,i) = 0.0
            UT_av(2,i) = 10.0*6.0*muL_av(i)/( U_av(1,i)*beta1*d(ncl)**2)
        case (-2)                            !-2 represents that the edge is the farfield  boundaries
            ux = U_av(2,i) - 0.5*( U_Rot(1,iedge(1,i)) + U_Rot(1,iedge(2,i)) )
            uy = U_av(3,i) - 0.5*( U_Rot(2,iedge(1,i)) + U_Rot(2,iedge(2,i)) )
            Vn = ux*vector(1,i)   + uy*vector(2,i)
            
            ! Vn = dot_product( U_av(2:3,i) , vector(:,i) ) 
            if ( Vn .LE. 0.0) then         !inflow
                UT_av(2,i) = C1*sqrt( u_inf**2 + v_inf **2 )/L
                UT_av(1,i) = muL_inf*10.0**(-C2)/rou_inf * UT_av(2,i)
            else
                UT_av(:,i) = UT(:,ncl)
            end if
            
        case default
        
           if( schemeType ==1 ) then
                if( Vn > 0.0)  then
                    UT_av(:,i)= UT(:,ncl)
                else
                    UT_av(:,i)= UT(:,ncr)  
                end if
            else
                UT_av(:,i)=0.5*( UT(:,ncl) + UT(:,ncr) )  
                
            end if
            
        end select
          
    end do
    
end subroutine