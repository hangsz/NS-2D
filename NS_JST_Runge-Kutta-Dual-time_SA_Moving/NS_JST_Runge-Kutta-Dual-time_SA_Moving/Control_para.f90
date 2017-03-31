module  Control_para    !Control parameters module
    implicit none
    
    !---------------------------------------------------
    !constant parameters
    real(8),parameter::pi=3.14159265
    !---------------------------------------------------
    !incoming flow information
    real(8),parameter::Ma_inf=0.6               !the mach number
    real(8),parameter::att=2.89                !attack angle,бу
    real(8),parameter::rou_inf=1.225            !density,kg/m^3
    real(8),parameter::p_inf=1.01325E5          !pressure,N/m^2
    real(8),parameter::muL_inf=1.78E-5        !molecular viscocity, pa s
    real(8),parameter::C=110.6                  !sutherland formula constant
    real(8),parameter::R=287.0                  !the ideal gas constant,J/(kg*K)      
    real(8),parameter::gamma=1.4                !the specific heat ratio
    real(8),parameter::PrL=0.72                 !Prantle coefficient
    real(8),parameter::PrT=0.9                 !Prantle coefficient
    !-------------------------------------------
    
    !Finite volume method,Jamson scheme
    integer,parameter::itermax=50          !the maximun iterative times
    real(8),parameter::k2=0.6,k4=1.0/64         !artifical dissipation,k2(1/2,1),k4(1/256,1/32)
    !real(8),parameter::CFL=2.0               !CFL number,explicit, less than 2*sqrt(2)
    real(8),parameter::CFL=1.0          !CFL number,implicit, 10E4~10E6
    real(8),parameter::eps=1.0E-10              !the convergence precision  of the average density
    
   !time step
    real(8),parameter::C_time =4.0
   
    !adjustable parameters of Runge_kutta scheme
    integer,parameter::stage=4                       !Runge_kutta's step number
    real(8),parameter::alpha(stage)=(/ 1.0/3, 4.0/15, 5.0/9, 1.0 /)    !optimised stage coefficients
    real(8),parameter::beta(stage)=(/ 1.0, 1.0/2, 0.0, 0.0 /)    !optimised stage coefficients
      !dual time
    integer,parameter::iter_inner =100
    !----------------------------------
    
     integer,parameter::phase=12           ! steps*100 per cycle,100
    real(8),parameter::dt_r = 0.00095      !real time step
    real(8),parameter::rot_cen(2)=(/0.273,0.0/)        !rotational center
    real(8),parameter::att_ampl = 2.41   !rotational amplitude,бу
    real(8),parameter::kr = 0.0808   !reduced frequency
    
    !----------------------------------

end module