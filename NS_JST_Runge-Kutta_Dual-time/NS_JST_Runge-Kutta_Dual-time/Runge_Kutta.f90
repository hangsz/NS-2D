subroutine Runge_Kutta
    implicit none
    integer::i
    integer::m,mm
    
  
    !****************************************************** 
    ! Multistage Schemes
    W0 = W  
    
    do m=1,stage
        call Mean_edge
        call Con_flux
        
        call Art_dissipation   
        
        call Gradient
        call Vis_flux
        
        if(m == 1)  then                 !calculate the dissipation and the time step only at the first stage 
            Rsi = Fc - Dissi - Fv 
            
        else
            Rsi =  Fc - (  beta(m)*Dissi + ( 1.0 - beta(m) ) * D_last  ) - Fv 
        end if 
       
        dt = CFL*vol/dt
        D_last = Dissi
        
        do mm= 1,5 
             Rsi(mm,:) = Rsi(mm,:) + 3.0/2/dt_r*vol*W(mm,:) - Q(mm,:)
             W(mm,:) = W0(mm,:) - alpha(m)*dt/vol  *  Rsi(mm,:) 
             
             !W(mm,:) = W0(mm,:) - alpha(m)*dt/vol * 1.0/(1.0 + 3.0/2/dt_r*alpha(m)*dt )  * ( Rsi(mm,:)- Q(mm,:) ) 
        end do
                
        U(1,:) = W(1,:)
        U(2,:) = W(2,:)/W(1,:)
        U(3,:) = W(3,:)/W(1,:)
        U(5,:) = (gamma-1)*( W(5,:)-U(1,:)*( U(2,:)**2+U(3,:)**2)/2.0 ) 
               
    end do
        
end subroutine



