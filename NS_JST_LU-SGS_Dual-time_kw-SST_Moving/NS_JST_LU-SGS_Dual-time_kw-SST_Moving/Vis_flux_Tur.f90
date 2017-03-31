subroutine   Vis_flux_Tur   !central difference
    implicit none
    
    integer::i
    integer::ncl,ncr                                  !the number of the left and the right cell
    real(8)::dUTdl(2)
    real(8)::UT_grad(2,2)
    real(8)::Fv_T(2)
    real(8)::sig_K,sig_omg
    !------------------------------------------------------------
    
     !set Fv to zero
     FvT=0.0            
     
     do i=1,nedges
        ncl=iedge(3,i)                      !atain the number of the left and the right cell
        ncr=iedge(4,i)
        
        !----------------------
        select case(ncr)
        case(-2:-1)
            
            dUTdl = (UT_av(:,i) - UT(:,ncl) ) / (0.5*lij(i) )  
            
            UT_grad = 1.0/2*Grad_UT(:,:,ncl)
        
            UT_grad(:,1) = UT_grad(:,1) - ( dot_product( UT_grad(:,1),tij(:,i) ) -dUTdl(1) ) * tij(:,i)
            UT_grad(:,2) = UT_grad(:,2) - ( dot_product( UT_grad(:,2),tij(:,i) ) -dUTdl(2) ) * tij(:,i)
             
            sig_K = 0.5*sigma_K(ncl) 
            sig_omg = 0.5*sigma_omg(ncl) 
        case default
        
            dUTdl = ( UT(:,ncr)-UT(:,ncl) ) / lij(i)
            
            UT_grad = 1.0/2*(Grad_UT(:,:,ncl) + Grad_UT(:,:,ncr) )
        
            UT_grad(:,1) = UT_grad(:,1)-( dot_product(UT_grad(:,1),tij(:,i) ) -dUTdl(1) ) * tij(:,i)
            UT_grad(:,2) = UT_grad(:,2)-( dot_product(UT_grad(:,2),tij(:,i) ) -dUTdl(2) ) * tij(:,i)
           
            sig_K = 0.5*(sigma_K(ncl)  + sigma_K(ncr) )
            sig_omg = 0.5*(sigma_omg(ncl)  + sigma_omg(ncr) )
        end select
        !-----------------------------------------------------
       
        Fv_T(1) = ( muL_av(i) + sig_K   * muT_av(i) )* ( vector(1,i) * UT_grad(1,1) + vector(2,i) * UT_Grad(2,1) )
        Fv_T(2) = ( muL_av(i) + sig_omg * muT_av(i) )* ( vector(1,i) * UT_grad(1,2) + vector(2,i) * UT_Grad(2,2) )
        
        if(ncr.GT. 0) then
            
            FvT(:,ncl) = FvT(:,ncl) + Fv_T
            FvT(:,ncr) = FvT(:,ncr) - Fv_T
        else
            FvT(:,ncl) = FvT(:,ncl) + Fv_T
        end if
        
     end do
  
end subroutine

    
    
    