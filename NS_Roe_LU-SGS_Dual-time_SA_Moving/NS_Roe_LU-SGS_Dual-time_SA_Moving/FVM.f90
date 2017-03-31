module FVM              
    !------------------------------------------------
        ! purpose:  the main body module to solve the N-S equations with the Finite Volume Method  
    !------------------------------------------------
    use controlPara
    use gridInfo
    use SA       
    
    implicit none   
    
    real(8)::u_inf,v_inf,a_inf                        
    real(8),allocatable,dimension(:,:):: U, W,           &
                                         U_av,           &
                                         Wn, Wn1, Q,     &        
                                               
                                         Fc,Fv,          &   
                                         ! Dissi,          &  
                                         Rsi                       

    real(8),allocatable,dimension(:)::   alf, alf_v,        &               
                                         lamda_c, lamda_v,  &        
                                         dt,                &               
                                         muL, muL_av
                                
    real(8),allocatable::Grad(:,:,:)
    
    !------------------------------
            !variables in relation to LU-SGS shceme
    integer,allocatable::assigned(:)
    integer,allocatable::lower(:,:),upper(:,:)            
    
    !------------------------------
          ! variables in relation to gird rigid moving
    real(8),allocatable::xy0(:,:),vector0(:,:),rij0(:,:),rR0(:,:),rL0(:,:),tij0(:,:)
    
    real(8),allocatable::U_Rot(:,:)
    
    real(8)::oldRsi=0.0
    
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
    ! oldRsi                     :      the mean value of the density-redisual
    !-----------------------------------------------------
    
contains

include "readFlow.f90"                     ! read stady flow information. 
include "rotation.f90"                     ! calculate the grid information afer rotating.
include "meanEdge.f90"                     ! calculate the flow properites on the edge.
include "conFlux_Roe.f90"                  ! calculate the convective flux.
include "timeStep.f90"                     !calcualate the time step.

!-------------------------------------------
  !Turbulence, get muT
   include "allocateMemory_Tur.f90"        ! allocate memory for the variables in relation to turbulence flow.
   include "meanEdge_tur.f90"              ! calculate the turbulent flow properties on the edge.
   include "conFlux_tur.f90"               ! calculate the convective flux of each cell in the turbulence control equation.
   include "artDissipation_tur.f90"        ! same as the above, artficial dissperation calcualtion.
   include "gradient_tur.f90"              ! the gradient of the turbulent flow properities
   include "visFlux_tur.f90"               ! viscous flux of each cell.
   include "SA_solver.f90"                 ! the main function in turbulence part.
   include "getmuT.f90"                    ! get the turbulent viscosity.
!-------------------------------------------
    
include "gradient.f90"                     ! the gradient of flow properities.  
include "visFlux.f90"                      ! the viscous flux of each cell.
               
include "solver.f90"                       ! the main function of FVM module.

include "reorder.f90"                      ! reorder the cell number sequence.
include "LU_SGS.f90"                       ! LU-SGS algorithm for algebraic  equations.

include "output.f90"                       ! output the flow solution when necessary.
include "outputFreeFlow.f90"               ! output the free flow information

subroutine allocateMemory
    !---------------------------------------------------
        ! purpose: allocate memory for all variables
    !---------------------------------------------------

    allocate( U(5,ncells) )        !the primative variables            
   
    allocate( W(5,ncells) )        !the conserved variables
    allocate( Wn(5,ncells) )      
    allocate( Wn1(5,ncells) )       
    allocate( Q(5,ncells) )        !the conserved variables
     
    allocate( U_av(6,nedges) )     !the primative variable on the edges, the boundary conditon is used here. 
    allocate( Grad(2,6,ncells) )   !the gradient of every cell,rou,u,v,w,p,T.
    allocate( Fc(5,ncells) )       !the convective flux
    ! allocate( Dissi(5,ncells) )   ! the artificial dissipation.
    
    ! viscosity
    allocate( Fv(5,ncells) )       !the viscous flux
    allocate( muL(ncells) )
    allocate( muL_av(nedges) )     
    ! time
    allocate( alf(nedges) )        !the spectrum radius of every edge
    allocate( alf_v(nedges) )
    allocate( lamda_c(ncells) )
    allocate( lamda_v(ncells))
    allocate( dt(ncells) )         ! the time  step
    
    allocate( Rsi(5,ncells) )      !the residual
    
    ! LU-SGS
    allocate( assigned(ncells) )   !the sequcence of the cells being reordered
    allocate( lower(3,ncells) )    !the lower edges of every cell,triangle has three edges
    allocate( upper(3,ncells) )    !the upper edges of every cell
    
    ! moving
    allocate( xy0(2,nnodes) )
    allocate( vector0(2,nedges) )
    allocate( rij0(2,nedges) ) 
    allocate( rR0(2,nedges) ) 
    allocate( rL0(2,nedges) )
    allocate( tij0(2,nedges) )
    
    allocate( U_Rot(2,nnodes) )
    
end subroutine

subroutine flowInit     
    !-------------------------------------------------
        ! purpose: initialize the flow field 
    !------------------------------------------------
    implicit none
   
    a_inf=sqrt(p_inf*gamma/rou_inf)
    u_inf=Ma_inf*a_inf*cosd(att)        
    v_inf=Ma_inf*a_inf*sind(att)
   

    W(1,:)=U(1,:)
    W(2,:)=U(1,:)*U(2,:)
    W(3,:)=U(1,:)*U(3,:)
    W(5,:)=U(5,:)/(gamma-1.0) + U(1,:)*(U(2,:)**2 + U(3,:)**2)/2.0  
    
    WT = U(1,:) * nuT
    
end subroutine


subroutine converge(flag)       
    !----------------------------------------------
        ! purpose:  verify wheather the flow solution converges
    !----------------------------------------------
    implicit none
    integer::i,j
    integer::flag                  !flag, 1:converge;0:disconverge
    real(8)::rou_ncell=0.0         !the mean density of n+1 layer
    real(8)::u_ncell=0.0           !the mean density of n+1 layer
    real(8)::v_ncell=0.0           !the mean density of n+1 layer
    real(8)::p_ncell=0.0           !the mean density of n+1 layer
      
    real(8),save::rou_mean = 1.225   !the mean density of n layer
    real(8),save::u_mean = 0.0       !the mean density of n layer
    real(8),save::v_mean = 0.0       !the mean density of n laye
    real(8),save::p_mean = 101325.0  !the mean density of n layer

    real(8)::newRsi=0.0

    !write(*,*)  "Converge"

    rou_ncell = sum(U(1,:))/ncells
    u_ncell = sum(U(2,:))/ncells
    v_ncell = sum(U(3,:))/ncells
    p_ncell = sum(U(5,:))/ncells
    
    flag = 0
    
    newRsi = sum(abs(Rsi(1,:)))/ncells

    if ( newRsi/oldRsi<1.0E-2)   flag = 1
    
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
    