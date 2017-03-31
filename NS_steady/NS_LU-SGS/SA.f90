module SA
    
    use Grid_info
    
    implicit none
    
    real(8),allocatable:: nu_SA(:),nu_av_SA(:)
    real(8),allocatable:: W_SA(:)
    
    real(8),allocatable:: Q_SA(:),dQdW_SA(:)
    
    real(8),allocatable:: Grad_nu_SA(:,:),Grad_W_SA(:,:),Fc_SA(:),Fv_SA(:),Dissi_SA(:),Rsi_SA(:)
    
    real(8),allocatable:: S_SA(:)
    real(8),allocatable:: chi_SA(:),fv1_SA(:),fv2_SA(:),fv3_SA(:)
    real(8),allocatable:: fw_SA(:),g_SA(:),rT_SA(:)
    
    real(8)::Cb1_SA = 0.1355,Cb2_SA = 0.6 , &
             Cv1_SA = 7.1,Cv2_SA = 5.0, sigma_SA =2.0/3, kapa_SA = 0.41 , &
             Cw1_SA , Cw2_SA =0.3 , Cw3_SA = 2.0 
             !Ct1 =1.0 ,Ct2 = 2.0,Ct3 = 1.3, Ct4 = 0.5   tansition
             
contains
    
!        turbulent model
!  turbulent model--SA
   include "Allocate_memory_SA.f90"
   include "Mean_edge_SA.f90"
   include "con_flux_SA.f90"
   include "Art_dissipation_SA.f90"
   include "Gradient_SA.f90"
   include "Vis_flux_SA.f90"
   include "Resi_SA.f90"
   include "getMut_SA.f90"
             
end module


    