subroutine kw_SST_solver
    implicit none
   
    integer::i
    real(8),allocatable::f1(:)
    
    real(8)::CDkw,arg1,arg2,P_T
    
    real(8)::txxF,tyyF,txyF
    
    allocate( f1(ncells) )
   
    !--------------------------------------------
    !write(*,*)   "Runge_kutta_Tur"
    
    !call Mean_edge_Tur
    call Gradient_Tur
      
    do i =1,ncells
        
        CDkw = max( 2.0*U(1,i) * sigma_omg2/UT(2,i) * ( Grad_UT(1,1,i)*Grad_UT(1,2,i) + Grad_UT(2,1,i)*Grad_UT(2,2,i) ) ,10.0**(-20)  )
        arg1 = min( max( sqrt(UT(1,i))/( 0.09*UT(2,i)*d(i) ) , 500.0*muL(i)/(U(1,i)*UT(2,i)*d(i)**2) ) , 4.0*U(1,i)*sigma_omg2 * UT(1,i)/ (CDkw*d(i)**2)  )
        
        f1(i) = tanh( arg1**4)
        
        arg2 = max( 2.0 *sqrt(UT(1,i))/( 0.09*UT(2,i)*d(i) ) , 500.0*muL(i)/(U(1,i)*UT(2,i)*d(i)**2) )
        f2(i) = tanh(arg2**2)
        
    end do
    !--------------------------
    !coefficients    
    sigma_K = f1*sigma_K1 + (1.0 -f1)*sigma_K2
    sigma_omg = f1*sigma_omg1 + (1.0 -f1)*sigma_omg2
    Cw = f1*Cw1 + (1.0 -f1)*Cw2
    beta = f1*beta1 + (1.0 -f1)*beta2
        
    !-------------------------
    call Con_flux_Tur  
    
    if( schemeType==2)  then
        call Art_dissipation_Tur
    else
        Dissi_T = 0
    end if
    
    call Vis_flux_Tur
    
    do i =1, ncells
        !txxF = 2.0*muT(i)*Grad(1,2,i) - ( 2.0*muT(i)/3.0) * ( Grad(1,2,i) + Grad(2,3,i) ) - 2.0/3*U(1,i)*UT(1,i)  
        !tyyF = 2.0*muT(i)*Grad(2,3,i) - ( 2.0*muT(i)/3.0) * ( Grad(1,2,i) + Grad(2,3,i) ) - 2.0/3*U(1,i)*UT(1,i)  
        !txyF = 2.0*muT(i)*0.5*( Grad(2,2,i)  + Grad(1,3,i) )
        !P_T = txxF*Grad(1,2,i) + tyyF*Grad(2,3,i) + txyF*( Grad(2,2,i) + Grad(1,3,i)  ) 
        !
        P_T = muT(i)*( Grad(2,2,i) - Grad(1,3,i))**2
        
        P_T = min(P_T,20.0*beta_star * U(1,i)*UT(2,i) *UT(1,i)  ) 
                                                                      
        QT(1,i) = P_T - beta_star * U(1,i)*UT(2,i) *UT(1,i)              
        QT(2,i) = Cw(i)*U(1,i)/muT(i)*P_T - beta(i)*U(1,i)*UT(2,i)**2 &
                &  + 2.0*(1.0 - f1(i))*U(1,i)*sigma_omg2/UT(2,i) *( Grad_UT(1,1,i)*Grad_UT(1,2,i) + Grad_UT(2,1,i)*Grad_UT(2,2,i) ) 
    end do
    
    !-----------------------------------------------
    do i=1,2    
        RsiT(i,:)= FcT(i,:) - Dissi_T(i,:) - FvT(i,:)  - QT(i,:)*vol
    end do
    
    
    dQTdWT(1,:) =  - beta_star  * UT(2,:)
    dQTdWT(2,:) =  - 2.0 * beta * UT(2,:)  
                    
    
end subroutine


                       