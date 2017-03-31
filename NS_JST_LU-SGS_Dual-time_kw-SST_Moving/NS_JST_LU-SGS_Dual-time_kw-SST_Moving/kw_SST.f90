module kw_SST
    implicit none
    
    real(8),allocatable:: muT(:),muT_av(:)
    
    real(8),allocatable:: UT(:,:) , UT_av(:,:)
    
    real(8),allocatable:: WT(:,:)
    
    real(8),allocatable:: WTn(:,:),WTn1(:,:),Q_WT(:,:)
    
    
    real(8),allocatable:: QT(:,:),dQTdWT(:,:)
    
    real(8),allocatable:: Grad_UT(:,:,:),Grad_WT(:,:,:),FcT(:,:),FvT(:,:),Dissi_T(:,:),RsiT(:,:)
    
    
    real(8)::C1 = 1.0,C2 =  5.0, L = 20.0 
    real(8)::a1 = 0.31,beta_star = 0.09 , kapa = 0.41
    
    
    real(8)::sigma_K1 = 0.85 , sigma_omg1 = 0.5 , beta1=0.075  , Cw1 = 0.533
    real(8)::sigma_K2 = 1.0  , sigma_omg2 = 0.856 , beta2=0.0828 , Cw2 = 0.440  
    
    real(8),allocatable::sigma_K(:), sigma_omg(:), beta(:), Cw(:) 
    real(8),allocatable::f2(:)
    
    !-------------------------------------------
    integer,parameter:: schemeType = 2      ! 1. upwind scheme
                                            ! 2. central difference scheme     
                                            
                                            
    
end module


    