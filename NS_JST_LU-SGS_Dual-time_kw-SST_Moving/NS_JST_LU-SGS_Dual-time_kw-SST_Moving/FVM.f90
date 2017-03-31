module FVM              !finite volume method 
    !object:  get the steady or unsteady flow properties finally 
   
    !use Control_para
    use Grid_info
    use kw_SST        !turbulence variables
    implicit none   
    real(8)::u_inf,v_inf,a_inf                        !incoming flow's velocities and sound velocity
    real(8),allocatable::U(:,:) ,&    !the oringinal variables of every cell
                         W(:,:)   , Wn(:,:) , Wn1(:,:),Q(:,:) ,&                 !the conservative variables of every cell
                         W0(:,:)  , &                 !the zeroth step conservative variables
                         Rsi(:,:) , &                 !residual 
                         Fc(:,:)  , &                 !convective flux of every cell
                         alf(:) ,alf_v(:), &                 !spectral radius of every edge
                         lamda_c(:) ,lamda_v(:), &   
                         Dissi(:,:)  , &              !artificial dissipation of every cell
                         dt(:)      ,&                !time step of every cell 
                         Grad(:,:,:)  ,&
                         Fv(:,:), muL(:), muL_av(:) ,&
                         U_av(:,:)          
    
    integer,allocatable::assigned(:)
    integer,allocatable::lower(:,:),upper(:,:)              !this two arrays are needed for LU-SGS method 
    
    
    !moving
    real(8),allocatable::xy0(:,:),vector0(:,:),rij0(:,:),rR0(:,:),rL0(:,:),tij0(:,:)
    
    real(8),allocatable::U_Rot(:,:)
    
contains

include "Read_flow.f90"
include "Rotation.f90"
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
   include "kw_SST_slover.f90"
   include "muT_cal.f90"
   
include "Gradient.f90"     
include "Vis_flux.f90"
  

include "Reorder.f90"
include "LU_SGS.f90"

include "Output.f90"
include "outputFreeFlow.f90"


subroutine Allocate_memory

    !allocate memory 
    allocate( U(5,ncells) )        !the primative variables            
   
    allocate( W(5,ncells) )        !the conserved variables
    allocate( Wn(5,ncells) )      
    allocate( Wn1(5,ncells) )       
    allocate( Q(5,ncells) )        !the conserved variables
     
    allocate( U_av(6,nedges) )     !the primative variable on the edges, the boundary conditon is used here 
    allocate( Grad(2,6,ncells) )   !the gradient of every cell,rou,u,v,w,p,T
    allocate( Fc(5,ncells) )       !the convective flux
    allocate( Dissi(5,ncells) )   ! the artificial dissipation
    
    !viscosity
    allocate( Fv(5,ncells) )       !the viscous flux
    allocate( muL(ncells) )
    allocate( muL_av(nedges) )     
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
    
     !moving
    allocate( xy0(2,nnodes) )
    allocate( vector0(2,nedges) )
    allocate( rij0(2,nedges) ) 
    allocate( rR0(2,nedges) ) 
    allocate( rL0(2,nedges) )
    allocate( tij0(2,nedges) )
    
    allocate( U_Rot(2,nnodes) )
    
end subroutine

subroutine Flow_init      !initialize the flow field

    implicit none
   
    a_inf=sqrt(p_inf*gamma/rou_inf)
    u_inf=Ma_inf*a_inf*cosd(att)        
    v_inf=Ma_inf*a_inf*sind(att)
    
    
    W(1,:)=U(1,:)
    W(2,:)=U(1,:)*U(2,:)
    W(3,:)=U(1,:)*U(3,:)
    W(5,:)=U(5,:)/(gamma-1.0) + U(1,:)*(U(2,:)**2 + U(3,:)**2)/2.0  + U(1,:)*UT(1,:)
    
    !initilize the turbulent flow
    WT(1,:) = U(1,:) * UT(1,:)
    WT(2,:) = U(1,:) * UT(2,:)
   
    
end subroutine

subroutine Solver          !the Solver
    implicit none
    integer::i,j
    integer::iter          !iterative variable
    integer::count
    integer::flag          !the variable to judge wheathe the mean density converges
    character(len=30)::filename
    real(8)::t_total = 0.0
    real(8)::omg,AoA,angular_rate
    
    write(*,*)  "Solver"
    
    call Grid
    call Allocate_memory  
    call Allocate_memory_Tur
    
    call outputFreeFlow
    call Read_flow
    call Flow_init
    
    call Reorder   !reordering is needed for LU-SGS
    
       
    !-------------------------------------------------
    !geometry set
    !translate the grid to where the rotational center's coordinates are zero
    
    vector0 = vector
    rij0 =rij
    rR0 = rR
    rL0 = rL
    tij0 = tij
    
    do i=1,nnodes
        xy0(:,i) = xy(:,i)-rot_cen(:)
    end do 
    !--------------------------------------
    !frequency
    omg =  2.0*kr*Ma_inf*a_inf
    
    !--------------------------------------
    
    Wn = W
    Wn1 = W
    
    WTn = WT
    WTn1 = WT
    
    do i=1,phase
    
        do iter = 1,itermax
             
            t_total = t_total + dt_r
            AoA = att + att_ampl* sin( omg * t_total )
            angular_rate = omg * ( att_ampl/180.0 *pi ) *cos( omg * t_total)
            
            Wn1 = Wn
            Wn = W
            
            WTn1 = WTn
            WTn = WT
                
            do j =1,5
                Q(j,:) = 2.0 /dt_r *vol*Wn(j,:) -1.0/2/dt_r * vol*Wn1(j,:)
            end do
            
            do j =1,2   
                Q_WT(j,:) = 2.0 /dt_r *vol*WTn(j,:) -1.0/2/dt_r * vol*WTn1(j,:)
            end do
            
            call Rotation(AoA,angular_rate)
         
             
            do count =1,iter_inner 
                write(*,*)  count , iter, i
                write(*,*) "t:",t_total,"dt:",dt_r
                write(*,*) "AoA:",AoA
          
                call Mean_edge       !calculate the primative variables on every edge
              
                !---------------------------------------------
                !the viscosity on the edges
                muL =   1.45*( U(5,:)/R/U(1,:) )**(3.0/2)/ ( U(5,:)/R/U(1,:) + 110.0 ) * 1.0E-6
                muL_av = 1.45*U_av(6,:)**(3.0/2)/( U_av(6,:)+110.0 ) *1.0E-6
              
                !turbulence
                call muT_cal
                !----------------------------------------------
                call Mean_edge_Tur   
                 
                call Con_flux        !calculate the convective flux using Roe scheme
        
                call Art_dissipation 
        
                dt = CFL*vol/(lamda_c + C_time*lamda_v )       !calculate the time step
        
                call Gradient        !calculate the gradient of every cell using the Mean_edge subroutine
           
                call kw_SST_slover       !calculate the muT
                
                do j=1,2
                    RsiT(j,:) = RsiT(j,:) + 3.0/2/dt_r *vol*WT(j,:) - Q_WT(j,:)
                end do
                
               
                
                call Vis_flux        !calculate the viscous flux
        
                Rsi = Fc - Dissi - Fv
             
                do j =1,5
                    Rsi(j,:) = Rsi(j,:) + 3.0/2/dt_r *vol*W(j,:) - Q(j,:)  
                end do
                
                call LU_SGS
        
                call Converge(flag)
        
                !write(*,*)
                !if(flag == 1)  then
                !    write(*,*)  "...Converge..."
                !    exit 
                !else
                !    write(*,*)  iter,"...Misconverge..."
                !end if    
            end do
        end do
        
    
        write(filename,"(I2)") i
            
        filename = "flow-info-" // trim(filename)
        call Output(filename) 
           
    end do
    
   
    
    
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
    