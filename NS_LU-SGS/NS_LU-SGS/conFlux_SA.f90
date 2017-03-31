subroutine conFlux_SA(U_av,U_Rot)        !convective flux
    implicit none
    real(8),intent(in)::U_av(:,:),U_Rot(:,:)     
    
    integer::i
    integer::ncl,ncr                                       !the number of the left and the right cell
                           
    real(8)::W_SA_av,ux,uy,Vn                                     
                              
    !set Fc_SA to zero
   
    Fc_SA = 0.0
    
    do i = 1,nedges
        ncl=iedge(3,i)
        ncr=iedge(4,i)
        
        ux =  U_av(2,i) - 0.5 * ( U_Rot(1,iedge(1,i)) + U_Rot(1,iedge(2,i)) )
        uy =  U_av(3,i) - 0.5 * ( U_Rot(2,iedge(1,i)) + U_Rot(2,iedge(2,i)) )
        Vn = ux* vector(1,i)   + uy * vector(2,i) 
      !  Vn= U_av(2,i)*vector(1,i) + U_av(3,i)*vector(2,i)
       
        W_SA_av = U_av(1,i)*nu_av_SA(i)
       
        if ( ncr .GT. 0) then
            !the left cell plus the convective flux
            
            Fc_SA(ncl)  = Fc_SA(ncl) + W_SA_av * Vn
            Fc_SA(ncr)  = Fc_SA(ncr) - W_SA_av * Vn
        else 
            Fc_SA(ncl)  = Fc_SA(ncl) + W_SA_av * Vn
        end if

    end do
    
end subroutine
