subroutine Resi_SA( U,U_av, alf,muL,muL_av ,Grad ,U_Rot)
    implicit none
   
    real(8),intent(in)::U(:,:),U_av(:,:),alf(:),muL(:),muL_av(:),Grad(:,:,:),U_Rot(:,:)
    
    real(8),allocatable::chi_max(:)
    
    allocate( chi_max(ncells) )
    
    !write(*,*)   "Resi_SA"
    
    !call Mean_edge_SA(U_av)
    call conFlux_SA( U_av,U_Rot ) 
    
    if ( schemeType==2)  then
        call artDissipation_SA( U, alf )
    else
        Dissi_SA = 0.0 
    end if
    
    call gradient_SA( U_av )
    call visFlux_SA( U_av,muL_av )
                            
   !---------------------------------------
    !coefficients 
      
    S_SA = abs( Grad(2,2,:) - Grad(1,3,:)  ) !omg11 =0,omg22=0,omg12**2= omg21**2
        
    chi_SA = nu_SA/( muL/U(1,:) ) 
    fv1_SA = chi_SA**3/( chi_SA**3 + Cv1_SA**3)
    fv2_SA = ( 1.0 + chi_SA/Cv2_SA)**(-3)
        
    !------------------------------
    where( chi_SA .GT. 0.001) 
        chi_max = chi_SA
    elsewhere 
        chi_max = 0.001
    end where
        
    fv3_SA = ( 1.0 + chi_SA*fv1_SA)*( 1.0 -fv2_SA)/ chi_max
    !----------------------------------------------
        
    S_SA = fv3_SA*S_SA + nu_SA/( kapa_SA*d)**2 * fv2_SA
        
        
    rT_SA = nu_SA / ( S_SA * kapa_SA**2 * d**2 )
        
    g_SA = rT_SA + Cw2_SA*( rT_SA**6 - rT_SA) 
        
    Cw1_SA=Cb1_SA/kapa_SA**2+(1+Cb2_SA)/sigma_SA
        
    fw_SA = g_SA* (1.0 + Cw3_SA**6)**(1.0/6) /( g_SA**6 + Cw3_SA**6 )**(1.0/6)
        
    Q_SA = Cb1_SA*S_SA*W_SA + Cb2_SA/sigma_SA*( Grad_W_SA(1,:) * Grad_nu_SA(1,:) + Grad_W_SA(2,:) * Grad_nu_SA(2,:) ) - Cw1_SA * fw_SA * U(1,:)*(nu_SA/d)**2
        
    !-----------------------------------------------------
      
    
    Rsi_SA= Fc_SA - Dissi_SA - Fv_SA  - Q_SA*vol
          

    !-----------------------------------------
    ! negative part is necessary to enhance the magnitude of the diaginal element.
    dQdW_SA =  Cb1_SA*S_SA  - 2.0*Cw1_SA*fw_SA*nu_SA/d**2
    
    where( dQdW_SA > 0.0)  
        dQdW_SA = 0.0
    end where
           
end subroutine


                       