subroutine SA_solver
    !------------------------------------------------
        ! purpose: the main function and start function of turbulent control equations 
    !------------------------------------------------
        
    implicit none
   
    real(8),allocatable::chi_max(:)
    real(8),allocatable::dQdWT(:)
    
    allocate( chi_max(ncells) )
    
    !write(*,*)   "Runge_kutta_Tur"
    
    call meanEdge_tur          ! get the turbulent varibales on the edge
    call conFlux_tur           ! calculate the convective flux
    call artDissipation_tur    ! calculate the artificial dissperation
    
    call gradient_tur          
    call visFlux_tur           ! calculate the viscous flux, the molecule viscosity is need
      
                         
    !coefficients 
      
    S = abs(Grad(2,2,:) - Grad(1,3,:) ) 
        
    chi = nuT/( muL/U(1,:) ) 
    fv1 = chi**3/( chi**3 + Cv1**3)
    fv2 = ( 1.0 + chi/Cv2)**(-3)
        
   
    where( chi .GT. 0.001) 
        chi_max = chi
    elsewhere 
        chi_max = 0.001
    end where
        
    fv3 = ( 1.0 + chi*fv1)*( 1.0 -fv2)/ chi_max
        
    S = fv3*S + nuT/( kapa*d)**2 * fv2
        
        
    rT = nuT / ( S * kapa**2 * d**2 )
        
    g = rT + Cw2*( rT**6 - rT) 
        
    Cw1=Cb1/kapa**2+(1+Cb2)/sigma
        
    fw = g* (1.0 + Cw3**6)**(1.0/6) /( g**6 + Cw3**6 )**(1.0/6)
        
    QT = Cb1*S*WT + Cb2/sigma*( Grad_WT(1,:) * Grad_nuT(1,:) + Grad_WT(2,:) * Grad_nuT(2,:) ) - Cw1 * fw  * U(1,:)*(nuT/d)**2
        
        
    RsiT= FcT - Dissi_T - FvT  - QT*vol
          
   
    dQTdWT =  Cb1*S  - 2.0*Cw1*fw*nuT/d**2
                    
    
end subroutine


                       