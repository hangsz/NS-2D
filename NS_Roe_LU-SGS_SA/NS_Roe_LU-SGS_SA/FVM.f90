module FVM            
    !------------------------------------------------
        ! purpose:  the main body module to solve the N-S equations with the Finite Volume Method  
    !------------------------------------------------
   
    !use Control_para
    use Grid_info
    use  SA                                         
    implicit none   
    real(8)::u_inf,v_inf,a_inf                    
    real(8),allocatable::U(:,:) ,&               
                         W(:,:)   , Wn(:,:) , Wn1(:,:),Q(:,:) ,&             
                         W0(:,:)  , &                
                         Rsi(:,:) , &          
                         Fc(:,:)  , &                
                         alf(:) ,alf_v(:), &         
                         lamda_c(:) ,lamda_v(:), &       
                         dt(:)      ,&             
                         Grad(:,:,:)  ,&
                         Fv(:,:), muL(:), muL_av(:) ,&
                         U_av(:,:)          
    
    integer,allocatable::assigned(:)
    integer,allocatable::lower(:,:),upper(:,:)            
    
    
    real(8),allocatable::U_Rot(:,:)
    
    !----------------------------------------------------
                    ! variable specification
                    
    ! u_inf,v_inf,a_inf          :      incoming flow's velocities and sound velocity, m/s
    ! U,W                        :      the primative variables and the conservative variables of the N-S equations of each cell.
    ! U_av                       :      the primative varibales of each edge
                                       ! U(1,:)  - density, kg/m^3
                                       ! U(2:3,:)  - the x-component and y-component of the velocity, m/s
                                       ! U(4,:)    - the alternative z-component of the velocity, m/s; it is no use in 2-D. 
                                       ! U(5,:)    - the static pressure, pa
                                       ! U(6,:)    - the temperature, K
                                       
    ! Wn,Wn1                     :      the intermediate conservative variables, it's used in dual time algorithm.
    ! Q                          :      the source term of dual time algorithm.
    
    ! Fc,Fv                      :      the convective and the viscous flux of the flow crossing the edge of each cell,using boundary flow properties. 
    ! Dissi                      :      the artificial dissipation of each cell.
    ! Rsi                        :      the residual
    
    ! alf,alf_v                  :      the convective and visous spectral radius of each edge.
    ! lamda_c,lamda_v            :      the total of alf and alf_v respectively of each cell.
    ! dt                         :      the local time step of each cell.
    ! muL,muL_av                 :      the molecule viscosity of each cell and each edge respectively.
    
    ! Grad                       :      the gradient of the variables of each cell.
    
    ! assigned                   :      stores the reorderd sequence of the cell number.
    ! lower,upper                :      stores the lower cells and upper cells of each cell.
    
    ! xy0,vector0,...            :      stores the initial grid information before rotating.         
    ! U_Rot                      :      rotating speed of each node.
    !-----------------------------------------------------
contains

include "Mean_edge.f90"
include "Gradient.f90" 

include "Con_flux_Roe.f90"
include "Time_step.f90"


 !Turbulence ,to get muT
   include "muT_cal.f90"
   include "Allocate_memory_Tur.f90"
   include "Mean_edge_Tur.f90" 
   include "Gradient_Tur.f90"
   
  ! include "con_flux_Roe_Tur.f90"
   include "Con_flux_Tur.f90"
   include "Art_dissipation_Tur.f90"
   
   include "Vis_flux_Tur.f90"
   include "SA_solver.f90"
   
   
include "Vis_flux.f90"
  
include "Reorder.f90"
include "LU_SGS.f90"

include "Output.f90"
include "Cp.f90"

subroutine Allocate_memory

    !allocate memory 
    allocate( U(5,ncells) )        !the primative variables            
   
    allocate( W(5,ncells) )        !the conserved variables
     
    allocate( U_av(6,nedges) )     !the primative variable on the edges, the boundary conditon is used here 
    allocate( Grad(2,6,ncells) )   !the gradient of every cell,rou,u,v,w,p,T
   
    allocate( alf(nedges) )
    allocate( alf_v(nedges))
    allocate( lamda_c(ncells) )
    allocate( lamda_v(ncells))
    allocate( Fc(5,ncells) )       !the convective flux
  
    !viscosity
    allocate( Fv(5,ncells) )       !the viscous flux
    allocate( muL(ncells) )
    allocate( muL_av(nedges) )     
    !time
    allocate( dt(ncells) )         ! the time  step
    
    allocate( Rsi(5,ncells) )      !the residual
    
    !LU-SGS
    allocate( assigned(ncells) )   !the sequcence of the cells being reordered
    allocate( lower(4,ncells) )    !the lower edges of every cell,triangle has three edges
    allocate( upper(4,ncells) )    !the upper edges of every cell
    
end subroutine

subroutine Flow_init      !initialize the flow field

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
    nuT = 0.1*muL_inf / rou_inf
    WT = rou_inf*nuT
    
end subroutine

subroutine Solver          !the Solver
    implicit none
    integer::i,j
    integer::iter          !iterative variable
    integer::count
    integer::flag          !the variable to judge wheathe the mean density converges
    character(len=30)::filename
    
    write(*,*)  "Solver"
    
    call Grid
    call Allocate_memory  
    call Allocate_memory_Tur
    call Flow_init
    
    call Reorder   !reordering is needed for LU-SGS
    
       
    do iter = 1,itermax
        write(*,*)  iter
        
        call Mean_edge       !calculate the primative variables on every edge
        !turbulence
        muL =   1.45*( U(5,:)/R/U(1,:) )**(3.0/2)/ ( U(5,:)/R/U(1,:) + 110.0 ) * 1.0E-6
        muL_av = 1.45*U_av(6,:)**(3.0/2)/( U_av(6,:)+110.0 ) *1.0E-6
       
        call muT_cal
        
        call Gradient        !calculate the gradient of every cell using the Mean_edge subroutine
       
        call Time_step
           
        call Con_flux        !calculate the convective flux using Roe scheme
      
        call SA_solver       !calculate the muT
        
        call Vis_flux        !calculate the viscous flux
        
        Rsi = Fc-Fv
        
        call LU_SGS
        
        call Converge(flag)
        
        if(flag == 1)  then
            write(*,*)  "...Converge..."
            exit 
        end if
        
    end do
    
    filename ="steady"
    call Output(filename)
    !write(*,*)  "Cp_cal"
    
end subroutine 


subroutine Converge(flag)          !verify wheather the flow converge
    implicit none
    integer::i,j
    integer::flag                  !flag, 1:converge;0:disconverge
    real(8)::rou_ncell=0.0    !the mean density of n+1 layer
    real(8)::u_ncell=0.0    !the mean density of n+1 layer
    real(8)::v_ncell=0.0    !the mean density of n+1 layer
    real(8)::p_ncell=0.0    !the mean density of n+1 layer
      
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
    write(*,*)  U(1,100),U(2,100),U(3,100),U(5,100) 
    write(*,*)  abs(rou_ncell-rou_mean),abs(p_ncell-p_mean)
    write(*,*)  abs(u_ncell-u_mean),abs(v_ncell-v_mean)
    write(*,*)  
      
    rou_mean = rou_ncell
    u_mean = u_ncell
    v_mean = v_ncell
    p_mean = p_ncell
    
end subroutine

end module
    