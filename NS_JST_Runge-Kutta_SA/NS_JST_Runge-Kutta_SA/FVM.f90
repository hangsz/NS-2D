module FVM              !finite volume method 
    !object:  get the steady or unsteady flow properties finally 
   
    use Control_para
    use Grid_info
    
    use  SA        !turbulence variables
    
    implicit none   
    real(8)::u_inf,v_inf,a_inf                        !incoming flow's velocities and sound velocity
    real(8),allocatable::U(:,:) ,U_av(:,:),&    !the oringinal variables of every cell
                         W(:,:)   , &                 !the conservative variables of every cell
                         W0(:,:)  , &                 !the zeroth step conservative variables
                         Rsi(:,:) , &                 !residual 
                         Fc(:,:)  , &                 !convective flux of every cell
                         alf(:)   , alf_v(:),&                 !spectral radius of every edge
                         Dissi(:,:)  , &              !artificial dissipation of every cell
                         dt(:)      ,&                !time step of every cell 
                         Grad(:,:,:)  ,&
                         Fv(:,:), muL(:), muL_av(:) 
            
                         
    !real(8)::dt_global                  
      
contains

include "Mean_edge.f90"
include "Con_flux.f90"
include "Art_dissipation.f90"

 !Turbulence ,to get muT
   include "Allocate_memory_Tur.f90"
   include "Mean_edge_Tur.f90"
   include "con_flux_Tur.f90"
   include "Art_dissipation_Tur.f90"
   include "Gradient_Tur.f90"
   include "Vis_flux_Tur.f90"
   include "SA_solver.f90"
   include "muT_cal.f90"
   
include "Gradient.f90"     
include "Vis_flux.f90"
  

include "Output.f90"
include "Cp.f90"

subroutine Allocate_memory

    !allocate memory 
    allocate( U(5,ncells) )        !the primative variables            
   
    allocate( W(5,ncells) )        !the conserved variables
    allocate( W0(5,ncells) )     !it's needed for Runge_kutta method
     
    allocate( U_av(6,nedges) )     !the primative variable on the edges, the boundary conditon is used here 
    allocate( Grad(2,6,ncells) )   !the gradient of every cell,rou,u,v,w,p,T
    allocate( Fc(5,ncells) )       !the convective flux
    allocate( Dissi(5,ncells) )   ! the artificial dissipation
    !time
    allocate( alf(nedges) )        !the spectrum radius of every edge
    allocate( alf_v(nedges) ) 
    allocate( dt(ncells) )         ! the time  step
    
    !viscosity
    allocate( Fv(5,ncells) )       !the viscous flux
    allocate( muL(ncells) )
    allocate( muL_av(nedges) )     
  
    allocate( Rsi(5,ncells) )      !the residual
    
   
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
   
    !initilize the turbulent viscosity
    nuT = 0.1* muL_inf / rou_inf
    WT=rou_inf*nuT
    
end subroutine

subroutine Solver          !the Solver
    implicit none
    integer::iter          !iterative variable
    integer::flag          !the variable to judge wheathe the mean density converges
    character(len=30)::filename
    
    write(*,*)  "Solver"
    
    call Grid
    call Allocate_memory   
    !turbulence
    call Allocate_memory_Tur
    call Flow_init
    
    do iter = 1,itermax
        write(*,*)  iter

        call Runge_kutta
        call Converge(flag)
        
        write(*,*) 
        if(flag == 1)  then
            write(*,*)  "...Converge..."
            exit 
        end if    
    end do
    
    write(*,*)  "Output"
    write(*,*)  "Cp_cal"
    filename ="steady"
    call Output(filename)
    !call Cp_cal(filename)
    
end subroutine 

subroutine Runge_kutta        !Runge-kutta scheme
    implicit none
    integer::m                !the step number of Runge-kutta
    integer::mm               !the number of equations
    character(len=30)::filename
    !write(*,*)   "Runge_kutta"
    !set w0 with the W
    W0 = W
    WT0 = WT
    do m = 1,Stage 
        
        call Mean_edge  
        muL =   1.45*( U(5,:)/R/U(1,:) ) **(3.0/2)/(U(5,:)/R/U(1,:) + 110.0) *1.0E-6
        muL_av = 1.45*U_av(6,:)**(3.0/2)/( U_av(6,:)+110.0 ) *1.0E-6
        call muT_cal
        
        call Con_flux  
        
        if(m .EQ. 1)  then                 !calculate the dissipation and the time step only at the first stage
            
            call Art_dissipation                     
            dt = CFL*vol/dt      !calculate the time step
         !   dt_global = minval( dt )
            
            call Gradient
            
            call SA_solver(m)       !calculate the muT

            call Vis_flux
    
        end if 
        
        do mm=1,5
            Rsi(mm,:)= Fc(mm,:) - Dissi(mm,:) - Fv(mm,:)
            W(mm,:)=W0(mm,:) - alpha(m)*dt * Rsi(mm,:) / vol    
        end do
      
        !calculate the original variables
        U(1,:) = W(1,:)
        U(2,:) = W(2,:)/W(1,:)
        U(3,:) = W(3,:)/W(1,:)
        U(5,:) = (gamma-1.0)*( W(5,:)-U(1,:)*( U(2,:)**2 + U(3,:)**2) / 2.0 )  
      
    end do 
    
    filename ="steady"
    call Output(filename)
    call Cp_cal(filename)
    
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
    real(8),save::p_mean = 103150.0  !the mean density of n layer


    !write(*,*)  "Converge"

    rou_ncell = sum(U(1,:))/ncells
    u_ncell = sum(U(2,:))/ncells
    v_ncell = sum(U(3,:))/ncells
    p_ncell = sum(U(5,:))/ncells
    
    flag = 0
    
    if (abs(rou_ncell-rou_mean) .LE. eps)   flag = 1
    
    write(*,*)  muT(1)
    write(*,*)  U(1,1),U(2,1),U(3,1),U(5,1)
    write(*,*)  abs(rou_ncell-rou_mean),abs(p_ncell-p_mean)
    write(*,*)  abs(u_ncell-u_mean),abs(v_ncell-v_mean)
    write(*,*)
       
    rou_mean = rou_ncell
    u_mean = u_ncell
    v_mean = v_ncell
    p_mean = p_ncell
    
end subroutine

end module
    