module  Control_para    !Control parameters module
    implicit none
    
    !------------------------------
    !constant parameters
    real(8),parameter::pi=3.14159265
    !-------------------------------
    
    !incoming flow information
    real(8),parameter::Ma_inf=0.6               !the mach number
    real(8),parameter::att=2.89                 !attack angle,бу
    real(8),parameter::rou_inf=1.225            !density,kg/m^3
    real(8),parameter::p_inf=1.01325E5          !pressure,N/m^2
    real(8),parameter::muL_inf=1.78E-5          !molecular viscocity, pa s
    real(8),parameter::C=110.6                  !sutherland formula constant
    real(8),parameter::R=287.0                  !the ideal gas constant,J/(kg*K)      
    real(8),parameter::gamma=1.4                !the specific heat ratio
    real(8),parameter::PrL=0.72                 !Prantle coefficient
    real(8),parameter::PrT=0.9                  !Prantle coefficient
    !------------------------------
     !time step
    real(8),parameter::C_time =4.0
    
    !Finite volume method,Jamson scheme
    integer,parameter::itermax=10000           !the maximun iterative times
    real(8),parameter::k2=0.6,k4=1.0/64         !artifical dissipation,k2(1/2,1),k4(1/256,1/32)
    !real(8),parameter::CFL=2.0                 !CFL num0.720ber,explicit, less than 2*sqrt(2)
    real(8),parameter::CFL=1000.0               !CFL number,implicit, 10E4~10E6
    real(8),parameter::eps=1.0E-10              !the convergence precision  of the average density
    
    !-----------------------------------------
    !limiter
     real(8)::K = 5.0
    !------------------------
    !-----------------------------
    !LU-SGS
    real(8)::omg = 1.5  !overrelaxation factor , 1< omg <= 2
    !-------------------------------------------
    character(len=1):: isNewGrid="y"                    ! determine wheather new gird is used.
    
    
    character(len=1),parameter:: gridNum = "1"   
    character(len=1),parameter:: schemeNum = "1"
    
    character(len=100)::caseRoute = "D:\Users\Hang\Documents\Visual Studio\Case"
    
end module