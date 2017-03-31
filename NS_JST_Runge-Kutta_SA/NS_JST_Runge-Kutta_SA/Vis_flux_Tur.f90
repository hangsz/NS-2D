subroutine   Vis_flux_Tur   !central difference
    implicit none
    
    integer::i
    integer::ncl,ncr                                  !the number of the left and the right cell
    real(8)::dnuTdl
    real(8)::nuT_grad(2)
    real(8)::txxT,tyyT
    real(8)::nuL_av
    !*************************************************************************
    
     !set Fv to zero
     FvT=0.0            
     
     do i=1,nedges
        ncl=iedge(3,i)                      !atain the number of the left and the right cell
        ncr=iedge(4,i)
        
        !************************************************************************
        select case(ncr)
        case(-2:-1)
            
            dnuTdl = (nuT_av(i) - nuT(ncl) ) / (0.5*lij(i) )  
            
            nuT_grad = 1.0/2*Grad_nuT(:,ncl)
        
            nuT_grad = nuT_grad - ( dot_product(nuT_grad,tij(:,i) ) -dnuTdl ) * tij(:,i)
                
        case default
        
            dnuTdl = ( nuT(ncr)-nuT(ncl) ) / lij(i)
            
            nuT_grad = 1.0/2*(Grad_nuT(:,ncl) + Grad_nuT(:,ncr) )
        
            nuT_grad = nuT_grad-( dot_product(nuT_grad,tij(:,i) ) -dnuTdl ) * tij(:,i)
            
        end select
        !----------------------------------------
        !tie molecure viscosity on the edge
        
        txxT = 1.0/sigma*( muL_av(i) + U_av(1,i)*nuT_av(i) )* nuT_grad(1)
        tyyT = 1.0/sigma*( muL_av(i) + U_av(1,i)*nuT_av(i) )* nuT_grad(2)
        
        if(ncr.GT. 0) then
            
            FvT(ncl) = FvT(ncl) + vector(1,i)*txxT + vector(2,i)*tyyT
            FvT(ncr) = FvT(ncr) -( vector(1,i)*txxT + vector(2,i)*tyyT )
        else
            FvT(ncl) = FvT(ncl) + vector(1,i)*txxT + vector(2,i)*tyyT
        end if
        
     end do
  
end subroutine

    
    
    