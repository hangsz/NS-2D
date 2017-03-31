subroUTine meanEdge_SST( u_inf,v_inf, U_av,muL_av ,U_Rot)  
    !object: get the flow properties on every edge, it's needed to calculate the gradient
    !        
    implicit none
    integer::i
    integer::ncl,ncr
    real(8),intent(in)::u_inf,v_inf
    real(8),intent(in)::U_av(:,:),muL_av(:),U_Rot(:,:)
    real(8)::ux,uy,Vn
    
    U_av_SST = 0.0
    
    do i=1,nedges
        ncl = iedge(3,i)                      ! the number of the left and the right cell
        ncr = iedge(4,i)
        select case( ncr ) 
        case (-1)                             !-1 represents that the edge is the aerofoil  surface         
            U_av_SST(1,i) = 0.0
            U_av_SST(2,i) = 10.0*6.0*muL_av(i)/( U_av(1,i)*beta1_SST*d(ncl)**2)
        case (-2)                            !-2 represents that the edge is the farfield  boundaries
          
            
            ux =  U_av(2,i) - 0.5 * ( U_Rot(1,iedge(1,i)) + U_Rot(1,iedge(2,i)) )
            uy =  U_av(3,i) - 0.5 * ( U_Rot(2,iedge(1,i)) + U_Rot(2,iedge(2,i)) )
            Vn = ux* vector(1,i)   + uy * vector(2,i)
            Vn = dot_product( U_av(2:3,i) , vector(:,i) ) 
          
            if ( Vn .LE. 0.0) then         !inflow
                U_av_SST(2,i) = C1_SST*sqrt( u_inf**2 + v_inf **2 )/L_SST
                U_av_SST(1,i) = muL_inf*10.0**(-C2_SST)/rou_inf * U_av_SST(2,i)
            else
                U_av_SST(:,i) = U_SST(:,ncl)
            end if
            
        case default
            if( schemeType ==1 ) then
                if( Vn > 0.0)  then
                    U_av_SST(:,i)= U_SST(:,ncl)
                else
                    U_av_SST(:,i)= U_SST(:,ncr)  
                end if
            else
                U_av_SST(:,i)=0.5*( U_SST(:,ncl) + U_SST(:,ncr) )  
                
            end if
            
        end select
          
    end do
    
end subroUTine
