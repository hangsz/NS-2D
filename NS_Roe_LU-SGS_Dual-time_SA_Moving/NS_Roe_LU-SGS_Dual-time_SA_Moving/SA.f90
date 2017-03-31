module SA
    !--------------------------------------------
        ! purpose: the variables and constants in relation to the SA turbulent model.
    !--------------------------------------------
    implicit none
    
    real(8),allocatable,dimension(:):: nuT, nuT_av,        &
                                       WT,                 &
                                       FcT, FvT,           &
                                       Dissi_T,            &
                                       QT, dQTdWT,         &
                                       RsiT,               &
    
                                       WTn, WTn1, Q_WT,    &     
                                       
                                       muT, muT_av,        &
                                       
                                       S,                  &
                                       chi, fv1, fv2, fv3, &
                                       fw, g, rT            
                                       
    real(8),allocatable,dimension(:,:):: Grad_nuT, Grad_WT
    
    
    real(8):: Cb1 = 0.1355, Cb2 = 0.6, &
              Cv1 = 7.1, Cv2 = 5.0, sigma =2.0/3, kapa = 0.41, &
              Cw1, Cw2 = 0.3, Cw3 = 2.0 
              ! Ct1 =1.0 ,Ct2 = 2.0,Ct3 = 1.3, Ct4 = 0.5   transition
             
    !--------------------------------------------------------------
                      ! variable specification
                      
    ! nuT, nuT_av               :       the primative variables in the cell center and on the edge.
    ! WT                        :       the conservative variable.
    ! FcT,FvT                   :       convective flux and viscous flux.
    ! Dissi_T                   :       artificial dissperation.
    ! QT,dQdWT                  :       the source term and the jacobian of dQ/dWT.
    ! RsiT                      :       the residual.
    
    ! WTn, WTn1, Q_WT           :       the variables are needed in dual time scheme.
    
    ! muT, muT_av               :       the turbulent viscosity in the cell center and on the edge.
    
    ! S                         :
    ! chi, fv1, fv2, fv3        :       variables relative to turbulent model.  
    ! fw, g, rT                 :
    
    !Grad_nuT, Grad_WT          :       the gradient of the primative and conservative variables in the cell center.
    
    !Cb1, Cb2, ...              :       constans in relation to turbulent model.
    !-------------------------------------------------------------- 
    
end module


    