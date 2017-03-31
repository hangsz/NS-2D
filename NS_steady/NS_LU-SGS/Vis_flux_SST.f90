subroUTine   Vis_flux_SST( U_av,muL_av,muT_av)  !central difference
    implicit none
    real(8),intent(in)::U_av(:,:),muL_av(:),muT_av(:)
    
    integer::i
    integer::ncl,ncr                                  !the number of the left and the right cell
    real(8)::dUTdl(2)
    real(8)::UT_grad(2,2)
    real(8)::Fv_T(2)
    real(8)::sig_K,sig_omg
    !------------------------------------------
    
     !set Fv to zero
     Fv_SST=0.0            
     
     do i=1,nedges
        ncl=iedge(3,i)                      !atain the number of the left and the right cell
        ncr=iedge(4,i)
        
        !-----------------------------------
        select case(ncr)
        case(-2:-1)
            
            dUTdl = (U_av_SST(:,i) - U_SST(:,ncl) ) / (0.5*lij(i) )  
            
            UT_grad = 1.0/2*Grad_UT(:,:,ncl)
        
            UT_grad(:,1) = UT_grad(:,1) - ( dot_product( UT_grad(:,1),tij(:,i) ) -dUTdl(1) ) * tij(:,i)
            UT_grad(:,2) = UT_grad(:,2) - ( dot_product( UT_grad(:,2),tij(:,i) ) -dUTdl(2) ) * tij(:,i)
            
            sig_K = sigma_K_SST(ncl)   
            sig_omg = sigma_omg_SST(ncl)   
            
        case default
        
            dUTdl = ( U_SST(:,ncr)-U_SST(:,ncl) ) / lij(i)
            
            UT_grad = 1.0/2*( Grad_UT(:,:,ncl) + Grad_UT(:,:,ncr) )
        
            UT_grad(:,1) = UT_grad(:,1)-( dot_product(UT_grad(:,1),tij(:,i) ) -dUTdl(1) ) * tij(:,i)
            UT_grad(:,2) = UT_grad(:,2)-( dot_product(UT_grad(:,2),tij(:,i) ) -dUTdl(2) ) * tij(:,i)
            
            sig_K   = 0.5*( sigma_K_SST(ncl)  +  sigma_K_SST(ncr)   )
            sig_omg = 0.5*( sigma_omg_SST(ncl)  +  sigma_omg_SST(ncr) )
        end select
        !--------------------------------------------------------
       
        
        Fv_T(1) = ( muL_av(i) + sig_K   * muT_av(i) )* ( vector(1,i) * UT_grad(1,1) + vector(2,i) * UT_grad(2,1) )
        Fv_T(2) = ( muL_av(i) + sig_omg * muT_av(i) )* ( vector(1,i) * UT_grad(1,2) + vector(2,i) * UT_grad(2,2) )
        
        if(ncr.GT. 0) then
            
            Fv_SST(:,ncl) = Fv_SST(:,ncl) + Fv_T
            Fv_SST(:,ncr) = Fv_SST(:,ncr) - Fv_T
        else
            Fv_SST(:,ncl) = Fv_SST(:,ncl) + Fv_T
        end if
        
     end do
  
end subroUTine

    
    
    