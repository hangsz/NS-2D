subroUTine   visFlux   !central difference
    implicit none
    
    integer::i
    integer::ncl,ncr                                  !the number of the left and the right cell
    real(8)::Ma                                       !mach number on any edge
    
    real(8)::T_av,T_inf  
               
    real(8)::c_av,mu_av,k_av,gradU(2),gradV(2),gradT(2)
    
    real(8)::lamda
    real(8)::txx,tyy,txy,theta_x,theta_y                             !txy,txx,tyy  every edge
    real(8)::dUdl,dVdl,dTdl
    
    !--------------------------------------------------
    
     !set Fv to zero
     Fv=0.0
     
     T_inf= p_inf/(R*rou_inf)
     
     call gradient                        !edge gradient is needed to calculate the viscous flux 
     
     do i=1,nedges
        ncl=iedge(3,i)                      !atain the number of the left and the right cell
        ncr=iedge(4,i)
        
        !--------------------------------------
        select case(ncr)
        case(-1)
            dUdl = ( U_av(2,i) - U(2,ncl) ) / ( 0.5 * lij(i) )
            dVdl = ( U_av(3,i) - U(3,ncl) ) / ( 0.5 * lij(i) )
            dTdl = ( U_av(6,i) - U(5,ncl)/U(1,ncl)/R ) / (0.5 * lij(i)  ) 
            
            gradU = dot_product( Grad(:,2,ncl),vector(:,i) ) *vector(:,i)/ds(i)**2         
            gradV = dot_product( Grad(:,3,ncl),vector(:,i) ) *vector(:,i)/ds(i)**2
            gradT = Grad(:,6,ncl) - dot_product( Grad(:,6,ncl),vector(:,i) )* vector(:,i)/ds(i)**2
        
            gradU = gradU-( dot_product(gradU,tij(:,i)) -dUdl ) * tij(:,i)
            gradV = gradV-( dot_product(gradV,tij(:,i)) -dVdl ) * tij(:,i)
            gradT = gradT-( dot_product(gradT,tij(:,i)) -dTdl ) * tij(:,i)
            
        case(-2)
            dUdl = ( U_av(2,i) - U(2,ncl) ) / ( 0.5 * lij(i) )
            dVdl = ( U_av(3,i) - U(3,ncl) ) / ( 0.5 * lij(i) )
            dTdl = ( U_av(6,i) - U(5,ncl)/U(1,ncl)/R ) / (0.5 * lij(i)  ) 
            
            gradU = 1.0/2*Grad(:,2,ncl)
            gradV = 1.0/2*Grad(:,3,ncl) 
            gradT = 1.0/2*Grad(:,6,ncl)
        
            gradU = gradU-( dot_product(gradU,tij(:,i)) -dUdl ) * tij(:,i)
            gradV = gradV-( dot_product(gradV,tij(:,i)) -dVdl ) * tij(:,i)
            gradT = gradT-( dot_product(gradT,tij(:,i)) -dTdl ) * tij(:,i)
     
        case default
            
            dUdl = ( U(2,ncr)-U(2,ncl) ) / lij(i)
            dVdl = ( U(3,ncr)-U(3,ncl) ) / lij(i)
            dTdl = ( U(5,ncr)/U(1,ncr)-U(5,ncl)/U(1,ncl) ) / R / lij(i) 
            
            gradU = 1.0/2*(Grad(:,2,ncl) + Grad(:,2,ncr) )
            gradV = 1.0/2*(Grad(:,3,ncl) + Grad(:,3,ncr) )
            gradT = 1.0/2*(Grad(:,6,ncl) + Grad(:,6,ncr) )
        
            gradU = gradU-( dot_product(gradU,tij(:,i)) -dUdl ) * tij(:,i)
            gradV = gradV-( dot_product(gradV,tij(:,i)) -dVdl ) * tij(:,i)
            gradT = gradT-( dot_product(gradT,tij(:,i)) -dTdl ) * tij(:,i)
            
        end select
        !---------------------------------------------------
       
        mu_av = muL_av(i) + muT_av(i)
        k_av= gamma*R/(gamma-1)*( muL_av(i)/PrL + muT_av(i)/PrT )       !  k=cp *mu/PrT  cp=gamma*R/gamma-1
        
        lamda = -2.0/3 *mu_av
        
        txx = lamda*( gradU(1)+gradV(2) ) + 2 * mu_av * gradU(1)
        tyy = lamda*( gradU(1)+gradV(2) ) + 2 * mu_av * gradV(2)
        txy = mu_av*( gradU(2)+gradV(1) )
        
        theta_x = U_av(2,i)*txx + U_av(3,i)*txy + k_av*gradT(1)
        theta_y = U_av(2,i)*txy + U_av(3,i)*tyy + k_av*gradT(2)
        
        
        Fv(2,ncl) = Fv(2,ncl) + vector(1,i)*txx + vector(2,i)*txy
        Fv(3,ncl) = Fv(3,ncl) + vector(1,i)*txy + vector(2,i)*tyy
        Fv(5,ncl) = Fv(5,ncl) + vector(1,i)*theta_x + vector(2,i)*theta_y
       
        if(ncr.GT. 0) then
            Fv(2,ncr) = Fv(2,ncr) - ( vector(1,i)*txx + vector(2,i)*txy )
            Fv(3,ncr) = Fv(3,ncr) - ( vector(1,i)*txy + vector(2,i)*tyy )
            Fv(5,ncr) = Fv(5,ncr) - ( vector(1,i)*theta_x + vector(2,i)*theta_y )
        end if
     end do
  
end subroUTine

    
    
    