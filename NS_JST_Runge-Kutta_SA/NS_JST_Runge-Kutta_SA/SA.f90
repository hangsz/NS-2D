module SA
    implicit none
    
    real(8),allocatable:: nuT(:),nuT_av(:)
    real(8),allocatable:: WT(:)  !rou*nuT
    real(8),allocatable:: WT0(:)
    
    real(8),allocatable:: muT(:),muT_av(:)
    
    real(8),allocatable:: QT(:)
    
    real(8),allocatable:: Grad_nuT(:,:),Grad_WT(:,:),FcT(:),FvT(:),Dissi_T(:),RsiT(:)
    
    real(8),allocatable:: S(:)
    real(8),allocatable:: chi(:),fv1(:),fv2(:),fv3(:)
    real(8),allocatable:: fw(:),g(:),rT(:)
    
    real(8)::Cb1 = 0.1355,Cb2 = 0.6 , &
             Cv1 = 7.1,Cv2 = 5.0, sigma =2.0/3, kapa = 0.41 , &
             Cw1 , Cw2 =0.3 , Cw3 = 2.0 
             !Ct1 =1.0 ,Ct2 = 2.0,Ct3 = 1.3, Ct4 = 0.5   tansition
end module


    