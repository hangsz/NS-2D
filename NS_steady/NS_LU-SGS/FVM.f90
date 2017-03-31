module FVM              !finite volume method 
    !object:  get the steady or unsteady flow properties finally 
   
    !use Grid_info
    use SA        
    use SST
    
    implicit none   
    real(8)::u_inf,v_inf,a_inf                        !incoming flow'S_SA velocities and sound velocity
    real(8),allocatable::U(:,:) , U_av(:,:) ,&                    !the oringinal variables of every cell
                         W(:,:)   , &                 !the conservative variables of every cell   
                         Fc(:,:)  , &                 !convective flux of every cell
                         alf(:) ,alf_v(:), &          !spectral radius of every edge
                         lamda_c(:) ,lamda_v(:), &    !spectral radius of every edge
                         Dissi(:,:)  , &              !artificial dissipation of every cell
                         dt(:)      ,&                !time step of every cell 
                         Grad(:,:,:)  ,&
                         Fv(:,:), &
                         Rsi(:,:) , &                 !residual 
                         muL(:), muL_av(:) 
                        
                        
                         
    real(8),allocatable:: muT(:),muT_av(:)
    
    integer,allocatable::assigned(:)
    
    integer,allocatable::lower(:,:),upper(:,:)              !this two arrays are needed for LU-SGS method 
    
contains

!---------------------------------------
! calculate the properties on the edge, considering the boundary conditions
include "Mean_edge.f90"
!----------------------------------------
       ! convective flux
!--------
! JST scheme
include "Con_flux_JST.f90"
include "Art_dissipation.f90"

!--------
! Roe scheme
include "Con_flux_Roe.f90"
include "Time_step.f90"
!-----------------------------------------
  
include "Gradient.f90"     
include "Vis_flux.f90"

!--------------------------------
! solver
include "solver_SA.f90"
include "solver_SST.f90" 
!---------------------------------
! algebra equation solver
include "Reorder.f90"
include "LU_SGS_SA.f90"
include "LU_SGS_SST.f90"

!-------------------------------
! post process
include "output_SA.f90"
include "output_SST.f90"
include "Cp.f90"
!-------------------------------

subroutine Allocate_memory

    !allocate memory 
    allocate( U(5,ncells) )        !the primative variables            
   
    allocate( W(5,ncells) )        !the conserved variables
     
    allocate( U_av(6,nedges) )     !the primative variable on the edges, the boundary conditon is used here 
    allocate( Grad(2,6,ncells) )   !the gradient of every cell,rou,u,v,w,p,T
    allocate( Fc(5,ncells) )       !the convective flux
    allocate( Dissi(5,ncells) )   ! the artificial dissipation
    
    !viscosity
    allocate( Fv(5,ncells) )       !the viscous flux
    allocate( muL(ncells) )
    allocate( muL_av(nedges) )     
    
    allocate( muT(ncells) ) 
    allocate( muT_av(nedges) )
    !time
    allocate( alf(nedges) )        !the spectrum radius of every edge
    allocate( alf_v(nedges) )
    allocate( lamda_c(ncells) )
    allocate( lamda_v(ncells))
    allocate( dt(ncells) )         ! the time  step
    
    allocate( Rsi(5,ncells) )      !the residual
    
    !LU-SGS
    allocate( assigned(ncells) )   !the sequcence of the cells being reordered
    allocate( lower(4,ncells) )    !the lower edges of every cell,triangle has three edges
    allocate( upper(4,ncells) )    !the upper edges of every cell
    
       
end subroutine

subroutine Flow_init_SA      !initialize the flow field

    implicit none
   
    a_inf=sqrt(p_inf*gamma/rou_inf)
    u_inf=Ma_inf*a_inf*cosd(att)        
    v_inf=Ma_inf*a_inf*sind(att)
   
    U(1,:)=rou_inf
    U(2,:)=u_inf
    U(3,:)=v_inf
    U(5,:)=p_inf
    
    W(1,:)=rou_inf
    W(2,:)=rou_inf*u_inf
    W(3,:)=rou_inf*v_inf
    W(5,:)=p_inf/(gamma-1.0) + rou_inf*(u_inf**2+v_inf**2)/2.0 
   
    !initilize the turbulent flow
    nu_SA = 0.1*muL_inf / rou_inf
    W_SA = rou_inf*nu_SA
    
   
end subroutine

subroutine Flow_init_SST      !initialize the flow field

    implicit none
   
    a_inf=sqrt(p_inf*gamma/rou_inf)
    u_inf=Ma_inf*a_inf*cosd(att)        
    v_inf=Ma_inf*a_inf*sind(att)
   

    U(1,:)=rou_inf
    U(2,:)=u_inf
    U(3,:)=v_inf
    U(5,:)=p_inf
    
    W(1,:)=rou_inf
    W(2,:)=rou_inf*u_inf
    W(3,:)=rou_inf*v_inf
    W(5,:)=p_inf/(gamma-1.0) + rou_inf*(u_inf**2+v_inf**2)/2.0    
    
    !initilize the turbulent flow

    muT = muL_inf*10.0**(-C2_SST)
    U_SST(2,:) = C1_SST*sqrt(u_inf**2 + v_inf**2)/ L_SST
    U_SST(1,:) = muL_inf*10.0**(-C2_SST)/rou_inf *U_SST(2,:)
    
    W_SST(1,:) = U(1,:) * U_SST(1,:)
    W_SST(2,:) = U(1,:) * U_SST(2,:)
    
    W(5,:) = W(5,:) + W_SST(1,:)   ! E= e + 0.5*v**2 + k
    
    
end subroutine

subroutine Converge(flag)          !verify wheather the flow converge
    implicit none
    integer::i,j
    integer::flag                  !flag, 1:converge;0:disconverge
    real(8)::rou_ncell=0.0         !the mean density of n+1 layer
    real(8)::u_ncell=0.0           !the mean density of n+1 layer
    real(8)::v_ncell=0.0           !the mean density of n+1 layer
    real(8)::p_ncell=0.0           !the mean density of n+1 layer
      
    real(8),save::rou_mean = 1.225  !the mean density of n layer
    real(8),save::u_mean = 0.0  !the mean density of n layer
    real(8),save::v_mean = 0.0  !the mean density of n laye
    real(8),save::p_mean = 101325.0  !the mean density of n layer


    !write(*,*)  "Converge"

    rou_ncell = sum(U(1,:))/ncells
    u_ncell = sum(U(2,:))/ncells
    v_ncell = sum(U(3,:))/ncells
    p_ncell = sum(U(5,:))/ncells
    
    flag = 0
    
    if (abs(rou_ncell-rou_mean) .LE. eps)   flag = 1
    write(*,*)  muT(100)
    write(*,*)  U(1,100),U(5,100)
    write(*,*)  U(2,100),U(3,100)
    write(*,*)  abs(rou_ncell-rou_mean),abs(p_ncell-p_mean)
    write(*,*)  abs(u_ncell-u_mean),abs(v_ncell-v_mean)
    write(*,*)  
      
    rou_mean = rou_ncell
    u_mean = u_ncell
    v_mean = v_ncell
    p_mean = p_ncell
    
end subroutine

end module
    