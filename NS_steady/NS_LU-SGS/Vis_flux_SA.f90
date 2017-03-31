subroutine   Vis_flux_SA( U_av,mul_av)   !central difference
    implicit none
    
    real(8),intent(in)::U_av(:,:),mul_av(:)
    integer::i
    integer::ncl,ncr                                  !the number of the left and the right cell
    real(8)::dnuTdl
    real(8)::nuT_grad(2)
    real(8)::txxT,tyyT
    real(8)::nuL_av
    !----------------------------------------------------------
    
     !set Fv to zero
     Fv_SA=0.0            
     
     do i=1,nedges
        ncl=iedge(3,i)                      !atain the number of the left and the right cell
        ncr=iedge(4,i)
        
        !------------------------------------------------
        select case(ncr)
        case(-2:-1)
            
            dnuTdl = (nu_av_SA(i) - nu_SA(ncl) ) / (0.5*lij(i) )  
            
            nuT_grad = 1.0/2*Grad_nu_SA(:,ncl)
        
            nuT_grad = nuT_grad - ( dot_product(nuT_grad,tij(:,i) ) -dnuTdl ) * tij(:,i)
                
        case default
        
            dnuTdl = ( nu_SA(ncr)-nu_SA(ncl) ) / lij(i)
            
            nuT_grad = 1.0/2*(Grad_nu_SA(:,ncl) + Grad_nu_SA(:,ncr) )
        
            nuT_grad = nuT_grad-( dot_product(nuT_grad,tij(:,i) ) -dnuTdl ) * tij(:,i)
            
        end select
        !-----------------------------------------------------------
        !tie molecure viscosity on the edge
        
        txxT = 1.0/sigma_SA*( muL_av(i) + U_av(1,i)*nu_av_SA(i) )* nuT_grad(1)
        tyyT = 1.0/sigma_SA*( muL_av(i) + U_av(1,i)*nu_av_SA(i) )* nuT_grad(2)
        
        if(ncr.GT. 0) then
            
            Fv_SA(ncl) = Fv_SA(ncl) + vector(1,i)*txxT + vector(2,i)*tyyT
            Fv_SA(ncr) = Fv_SA(ncr) -( vector(1,i)*txxT + vector(2,i)*tyyT )
        else
            Fv_SA(ncl) = Fv_SA(ncl) + vector(1,i)*txxT + vector(2,i)*tyyT
        end if
        
     end do
  
end subroutine

    
    
    