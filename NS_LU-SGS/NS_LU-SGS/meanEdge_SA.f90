subroutine meanEdge_SA( U_av,U_Rot)  
    !object: get the flow properties on every edge, it'S_SA needed to calculate the gradient
    !        
    implicit none
    integer::i
    integer::ncl,ncr
    real(8),intent(in)::U_av(:,:),U_Rot(:,:)
    real(8)::ux,uy,Vn
    
    nu_av_SA = 0.0
    
    do i=1,nedges
        ncl = iedge(3,i)                      ! the number of the left and the right cell
        ncr = iedge(4,i)
        select case( ncr ) 
        case (-1)                             !-1 represents that the edge is the aerofoil  surface         
            nu_av_SA(i) = 0.0
        case (-2)  !-2 represents that the edge is the farfield  boundaries
            
            ux =  U_av(2,i) - 0.5 * ( U_Rot(1,iedge(1,i)) + U_Rot(1,iedge(2,i)) )
            uy =  U_av(3,i) - 0.5 * ( U_Rot(2,iedge(1,i)) + U_Rot(2,iedge(2,i)) )
            Vn = ux* vector(1,i)   + uy * vector(2,i)
            !Vn = dot_product( U_av(2:3,i) , vector(:,i) ) 
            if ( Vn .LE. 0.0) then                     
                nu_av_SA(i)= 0.1 * muL_inf / rou_inf
            else
                nu_av_SA(i) = nu_SA(ncl)
            end if
            
        case default
        
           if( schemeType ==1)  then
           !  one order upwind scheme. the reasons to use it are that the convective effect is minor and the  scheme is more stable. 
               Vn = dot_product( U_av(2:3,i) , vector(:,i) ) 
               if( Vn > 0.0)  then              ! the direction of the vector is clockwise
                   nu_av_SA(i) = nu_SA(ncl)
               else
                   nu_av_SA(i) = nu_SA(ncr)
               end if
           else
           ! central difference scheme
               nu_av_SA(i)=0.5*( nu_SA(ncl) + nu_SA(ncr) )  
           end if
           
        end select
          
    end do
    
end subroutine