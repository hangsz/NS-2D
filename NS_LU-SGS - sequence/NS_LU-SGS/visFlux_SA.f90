subroutine   visFlux_SA( U_av,mul_av)   !central difference
    implicit none
    
    real(8),intent(in)::U_av(:,:),mul_av(:)
    integer::i
    integer::ncl,ncr                                  !the number of the left and the right cell
    real(8)::dnuTdl
    real(8)::nu_grad(2)
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
        case(-1)
            
            dnuTdl = (U_av_SA(i) - U_SA(ncl) ) / (0.5*lij(i) )  
            
            nu_grad = dot_product( Grad_U_SA(:,ncl),vector(:,i) ) * vector(:,i)/ds(i)**2
        
            nu_grad = nu_grad - ( dot_product(nu_grad,tij(:,i) ) -dnuTdl ) * tij(:,i)
            
        case(-2)
            
            dnuTdl = (U_av_SA(i) - U_SA(ncl) ) / (0.5*lij(i) )  
            
            nu_grad = 1.0/2*Grad_U_SA(:,ncl)
        
            nu_grad = nu_grad - ( dot_product(nu_grad,tij(:,i) ) -dnuTdl ) * tij(:,i)
                
        case default
        
            dnuTdl = ( U_SA(ncr)-U_SA(ncl) ) / lij(i)
            
            nu_grad = 1.0/2*(Grad_U_SA(:,ncl) + Grad_U_SA(:,ncr) )
        
            nu_grad = nu_grad-( dot_product(nu_grad,tij(:,i) ) -dnuTdl ) * tij(:,i)
            
        end select
        !-----------------------------------------------------------
        !tie molecure viscosity on the edge
        
        txxT = 1.0/sigma_SA*( muL_av(i) + U_av(1,i)*U_av_SA(i) )* nu_grad(1)
        tyyT = 1.0/sigma_SA*( muL_av(i) + U_av(1,i)*U_av_SA(i) )* nu_grad(2)
        
        Fv_SA(ncl) = Fv_SA(ncl) +  vector(1,i)*txxT + vector(2,i)*tyyT
        if(ncr.GT. 0) then
            Fv_SA(ncr) = Fv_SA(ncr) - ( vector(1,i)*txxT + vector(2,i)*tyyT )
        end if
        
     end do
  
end subroutine

    
    
    