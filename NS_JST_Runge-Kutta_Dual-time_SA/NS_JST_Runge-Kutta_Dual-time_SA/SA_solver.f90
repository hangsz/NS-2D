subroutine SA_solver(m)
   
    implicit none
    integer::m                !the step number of Runge-kutta
    integer::mm               !the number of equations
    
   
    real(8),allocatable::chi_max(:)
    real(8),allocatable::dQdWT(:)
    
   
    allocate( chi_max(ncells) )
    allocate( dQdWT(ncells) ) 
    
    !write(*,*)   "Runge_kutta_Tur"
    
    call Mean_edge_Tur
    call Con_flux_Tur 
    call Art_dissipation_Tur               
    call Gradient_Tur
    call Vis_flux_Tur
      
                    !calculate the dissipation and the time step only at the first stage
                                
    !************************
    !coefficients 
      
    S = abs(Grad(2,2,:) - Grad(1,3,:) ) !omg11 =0,omg22=0,omg12**2= omg21**2
        
    chi = nuT/( muL/U(1,:) ) 
    fv1 = chi**3/( chi**3 + Cv1**3)
    fv2 = ( 1.0 + chi/Cv2)**(-3)
        
    !***
    where( chi .GT. 0.001) 
        chi_max = chi
    elsewhere 
        chi_max = 0.001
    end where
        
    fv3 = ( 1.0 + chi*fv1)*( 1.0 -fv2)/ chi_max
    !***
        
    S = fv3*S + nuT/( kapa*d)**2 * fv2
        
        
    rT = nuT / ( S * kapa**2 * d**2 )
        
    g = rT + Cw2*( rT**6 - rT) 
        
    Cw1=Cb1/kapa**2+(1+Cb2)/sigma
        
    fw = g* (1.0 + Cw3**6)**(1.0/6) /( g**6 + Cw3**6 )**(1.0/6)
        
    QT = Cb1*S*WT + Cb2/sigma*( Grad_WT(1,:) * Grad_nuT(1,:) + Grad_WT(2,:) * Grad_nuT(2,:) ) - Cw1 * fw * U(1,:)* (nuT/d)**2
        
        
    !****************************
    if(m == 1)  then                 !calculate the dissipation and the time step only at the first stage 
            RsiT= FcT - Dissi_T - FvT  -QT*vol
            
    else
        RsiT =  FcT - (  beta(m)*Dissi_T + ( 1.0 - beta(m) ) * D_last_T  ) - FvT  -QT*vol
    end if 
       
     
    D_last_T = Dissi_T
        
    RsiT = RsiT + 3.0/2/dt_r*vol*WT - Q_WT
      
    WT = WT0 - alpha(m)*dt/vol  *  RsiT
    
    !W(mm,:) = W0(mm,:) - alpha(m)*dt/vol * 1.0/(1.0 + 3.0/2/dt_r*alpha(m)*dt )  * ( Rsi(mm,:)- Q(mm,:) )   
    ! dQdWT =  Cb1*S  - 2.0*Cw1*fw*nuT/d**2   
    ! WT= WT0 - alpha(m)/vol * 1.0/( 1.0/dt_global - dQdWT ) * RsiT 
        
    nuT = WT/U(1,:)
        
    !if the nu less than zero,set it to a small positive number
    where( nuT .LT.0.0 ) 
        nuT = minval( abs(nuT))
    end where
         
     
end subroutine

                       