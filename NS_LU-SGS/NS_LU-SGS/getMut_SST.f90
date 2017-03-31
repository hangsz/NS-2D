subroUTine  getMUT_SST( U,U_av, muL, muT, muT_av,U_Rot)
    implicit none
    real(8),intent(in)::U(:,:),U_av(:,:),muL(:),muT(:),U_Rot(:,:)
    real(8),intent(out)::muT_av(:)
    
    integer::i
    integer::iter
    integer::flag
    integer::ncl,ncr
    real(8)::ux,uy,Vn
    
    do i=1,nedges
        ncl = iedge(3,i)                      !atain the number of the left and the right cell
        ncr = iedge(4,i)
        select case( ncr ) 
        case (-1)                           !-1 represents that the edge is the aerofoil  surface         
            muT_av(i) = 0.0
        case (-2)   !-2 represents that the edge is the farfield  boundaries
            
            ux =  U_av(2,i) - 0.5 * ( U_Rot(1,iedge(1,i)) + U_Rot(1,iedge(2,i)) )
            uy =  U_av(3,i) - 0.5 * ( U_Rot(2,iedge(1,i)) + U_Rot(2,iedge(2,i)) )
            Vn = ux* vector(1,i)   + uy * vector(2,i)
            !Vn = dot_product( U_av(2:3,i) , vector(:,i) ) 
            if ( Vn .LE. 0.0) then
                muT_av(i)= muL_inf*10.0**(-C2_SST)
            else
                muT_av(i) = muT(ncl)
            end if
            
        case default
        
             muT_av(i)=0.5*( muT(ncl) + muT(ncr) )  
            
        end select
     
    end do
    
end subroUTine
