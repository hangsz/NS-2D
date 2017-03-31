subroUTine Resi_SST( u_inf,v_inf,U,U_av,muL,muL_av,muT,muT_av,Grad,U_Rot)
    implicit none
    real(8),intent(in)::u_inf,v_inf,U(:,:),U_av(:,:),muL(:),muL_av(:),muT(:),muT_av(:),Grad(:,:,:),U_Rot(:,:)
    
    integer::i
    real(8),allocatable::f1(:)
    
    real(8)::CDkw,arg1,arg2,P_T
    
    real(8)::txxF,tyyF,txyF
    
    allocate( f1(ncells) )
   
    !--------------------------------------------
    !write(*,*)   "Runge_kUTta_Tur"
    
    call meanEdge_SST( u_inf,v_inf,U_av,muL_av,U_Rot) 
    call gradient_SST
      
    do i =1,ncells
        
        CDkw = max( 2.0*U(1,i) * sigma_omg2_SST/U_SST(2,i) * ( Grad_SST(1,1,i)*Grad_SST(1,2,i) + Grad_SST(2,1,i)*Grad_SST(2,2,i) ) ,10.0**(-20)  )
        arg1 = min( max( sqrt(U_SST(1,i))/( 0.09*U_SST(2,i)*d(i) ) , 500.0*muL(i)/(U(1,i)*U_SST(2,i)*d(i)**2) ) , 4.0*U(1,i)*sigma_omg2_SST * U_SST(1,i)/ (CDkw*d(i)**2)  )
        
        f1(i) = tanh( arg1**4)
        
        arg2 = max( 2.0 *sqrt(U_SST(1,i))/( 0.09*U_SST(2,i)*d(i) ) , 500.0*muL(i)/(U(1,i)*U_SST(2,i)*d(i)**2) )
        f2_SST(i) = tanh(arg2**2)
        
    end do
    !--------------------------
    !coefficients    
    sigma_K_SST = f1*sigma_K1_SST + (1.0 -f1)*sigma_K2_SST
    sigma_omg_SST = f1*sigma_omg1_SST + (1.0 -f1)*sigma_omg2_SST
    Cw_SST = f1*Cw1_SST + (1.0 -f1)*Cw2_SST
    beta_SST = f1*beta1_SST + (1.0 -f1)*beta2_SST
        
    !-------------------------
    call conFlux_SST( U_av,U_Rot) 
    
    if( schemeType==2)  then
        call artDissipation_SST(U)
    else
        Dissi_SST = 0
    end if
    
    call visFlux_SST(muL_av,muT_av)
    
    do i =1, ncells
        !txxF = 2.0*muT(i)*Grad(1,2,i) - ( 2.0*muT(i)/3.0) * ( Grad(1,2,i) + Grad(2,3,i) ) - 2.0/3*U(1,i)*U_SST(1,i)  
        !tyyF = 2.0*muT(i)*Grad(2,3,i) - ( 2.0*muT(i)/3.0) * ( Grad(1,2,i) + Grad(2,3,i) ) - 2.0/3*U(1,i)*U_SST(1,i)  
        !txyF = 2.0*muT(i)*0.5*( Grad(2,2,i)  + Grad(1,3,i) )
        !P_T = txxF*Grad(1,2,i) + tyyF*Grad(2,3,i) + txyF*( Grad(2,2,i) + Grad(1,3,i)  ) 
        ! simplify the P_T
        P_T = muT(i)*( Grad(2,2,i) - Grad(1,3,i))**2
        
        P_T = min(P_T,20.0*beta_star_SST * U(1,i)*U_SST(2,i) *U_SST(1,i)  ) 
                                                                      
        Q_SST(1,i) = P_T - beta_star_SST * U(1,i)*U_SST(2,i) *U_SST(1,i)              
        Q_SST(2,i) = Cw_SST(i)*U(1,i)/muT(i)*P_T - beta_SST(i)*U(1,i)*U_SST(2,i)**2 &
                &  + 2.0*(1.0 - f1(i))*U(1,i)*sigma_omg2_SST/U_SST(2,i) *( Grad_SST(1,1,i)*Grad_SST(1,2,i) + Grad_SST(2,1,i)*Grad_SST(2,2,i) ) 
    end do
    
    !-----------------------------------------------
    do i=1,2    
        Rsi_SST(i,:)= Fc_SST(i,:) - Dissi_SST(i,:) - Fv_SST(i,:)  - Q_SST(i,:)*vol
    end do
    
    
    dQdW_SST(1,:) =  - beta_star_SST  * U_SST(2,:)
    dQdW_SST(2,:) =  - 2.0 * beta_SST * U_SST(2,:)  
                    
    
end subroUTine


                       