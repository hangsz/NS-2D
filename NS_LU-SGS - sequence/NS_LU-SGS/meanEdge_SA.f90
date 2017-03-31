subroutine meanEdge_SA( U_av,U_Rot)  
    !object: get the flow properties on every edge, it'S_SA needed to calculate the gradient
    !        
    implicit none
    integer::i
    integer::ncl,ncr
    real(8),intent(in)::U_av(:,:),U_Rot(:,:)
    real(8)::ux,uy,Vn
    
    U_av_SA = 0.0
    
    do i=1,nedges
        ncl = iedge(3,i)                      ! the number of the left and the right cell
        ncr = iedge(4,i) 
        
        ux =  U_av(2,i) - 0.5 * ( U_Rot(1,iedge(1,i)) + U_Rot(1,iedge(2,i)) )
        uy =  U_av(3,i) - 0.5 * ( U_Rot(2,iedge(1,i)) + U_Rot(2,iedge(2,i)) )
        Vn = ux* vector(1,i)/ds(i)   + uy * vector(2,i)/ds(i)
        !Vn = dot_product( U_av(2:3,i) , vector(:,i) ) 
        select case( ncr ) 
        case (-1)                             !-1 represents that the edge is the aerofoil  surface         
            U_av_SA(i) = 0.0
        case (-2)  !-2 represents that the edge is the farfield  boundaries
            if ( Vn .LE. 0.0) then                     
                U_av_SA(i)= 0.1 * muL_inf / rou_inf
            else
                U_av_SA(i) = U_SA(ncl)
            end if
            
        case default
        
            if( schemeType ==1)  then
            !  one order upwind scheme. the reasons to use it are that the convective effect is minor and the  scheme is more stable. 
                Vn = dot_product( U_av(2:3,i) , vector(:,i) ) 
                if( Vn > 0.0)  then              ! the direction of the vector is clockwise
                    U_av_SA(i) = U_SA(ncl)
                else
                    U_av_SA(i) = U_SA(ncr)
                end if
            else
            ! central difference scheme
                U_av_SA(i)=0.5*( U_SA(ncl) + U_SA(ncr) )  
            end if
           
        end select
          
    end do
    
end subroutine