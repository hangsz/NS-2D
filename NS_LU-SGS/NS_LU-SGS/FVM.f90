module FVM              !finite volume method 
    ! purpose:  get the steady or unsteady flow properties finally 
   
    !use gridInfo
    use SA        
    use SST
    use moving
    
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
                        
    !-----------------
    ! dual time
    real(8)::dt_r
    real(8),allocatable::Wn(:,:),Wn1(:,:),QW(:,:) 
    real(8)::oldRsi = 0.0    
    !--------------
    real(8),allocatable:: muT(:),muT_av(:)
    
    integer,allocatable::assigned(:)
    
    integer,allocatable::lower(:,:),upper(:,:)              !this two arrays are needed for LU-SGS method 
    
    !--------------------------
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
    ! QW                         :      the source term of dual time algorithm.
    
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
    
    ! U_Rot                      :      rotating speed of each node.
    
    !-----------------------------------------------------
contains
!----------------------------------
include "flowInit.f90"
include "flowInit_static.f90"
include "flowInit_dynamic.f90"
!---------------------------------------
! calculate the properties on the edge, considering the boundary conditions
include "meanEdge.f90"
!----------------------------------------
       ! convective flux
!--------
! JST scheme
include "conFlux_JST.f90"
include "artDissipation.f90"

!--------
! Roe scheme
include "conFlux_Roe.f90"
include "timeStep.f90"
!-----------------------------------------
       ! viscous flux
include "gradient.f90"     
include "visFlux.f90"

!--------------------------------
! solver
include "solver.f90"
include "solver_static.f90" 
include "solver_dynamic.f90" 
!---------------------------------
! algebra equation solver
include "reorder.f90"
include "LU_SGS.f90"
include "LU_SGS_SA.f90"
include "LU_SGS_SST.f90"

!-------------------------------
! post process
include "outputFreeFlow.f90"
include "output.f90"

include "Cp.f90"
!-------------------------------

subroutine allocateMemory

    write(*,*)  "allocateMemory"

    !allocate memory 
    allocate( U(5,ncells) )        !the primative variables            
   
    allocate( W(5,ncells) )        !the conserved variables
    !---------------------------
    ! dual time
    allocate( U_Rot(2,nnodes ) )
    if( moveNum /= 0 ) then
        
        allocate( Wn(5,ncells) )        !the conserved variables
        allocate( Wn1(5,ncells) )        !the conserved variables
        allocate( QW(5,ncells) )        !the conserved variables
    end if
    
    !-------------------------
    
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
    real(8),save::newRsi =0.0

    !write(*,*)  "Converge"

    rou_ncell = sum(U(1,:))/ncells
    u_ncell = sum(U(2,:))/ncells
    v_ncell = sum(U(3,:))/ncells
    p_ncell = sum(U(5,:))/ncells
    
    flag = 0
    !if( moveNum==0)  then
    !    if (abs(rou_ncell-rou_mean) .LE. eps)   flag = 1
    !else
        newRsi = sum(abs(Rsi(1,:)))/ncells
        if (newRsi/oldRsi < eps) flag = 1
    !end if
    
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
    