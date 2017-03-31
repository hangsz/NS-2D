module  Control_para    !Control parameters module
    implicit none
    
    !------------------------------
    !constant parameters
    real(8),parameter::pi=3.14159265
    !-------------------------------
    
    !incoming flow information
    real(8),parameter::Ma_inf=0.6               !the mach number
    real(8),parameter::att=2.89                 !attack ang_SAle,бу
    real(8),parameter::rou_inf=1.225            !density,kg_SA/m^3
    real(8),parameter::p_inf=1.01325E5          !pressure,N/m^2
    real(8),parameter::muL_inf=1.78E-5          !molecular viscocity, pa S_SA
    real(8),parameter::C=110.6                  !sutherland formula constant
    real(8),parameter::R=287.0                  !the ideal gas constant,J/(kg_SA*K)      
    real(8),parameter::gamma=1.4                !the specific heat ratio
    real(8),parameter::PrL=0.72                 !Prantle coefficient
    real(8),parameter::PrT=0.9                  !Prantle coefficient
    !------------------------------
    !time step
    real(8),parameter::C_time =4.0
    
    !------------------------------------
    !Finite volume method , Jamson scheme
    integer,parameter::itermax=10000            !the maximun iterative times
    real(8),parameter::k2=0.6,k4=1.0/64         !arT_SAifical dissipation,k2(1/2,1),k4(1/256,1/32)
    !real(8),parameter::CFL=2.0                 !CFL number,explicit, less than 2*sqrt(2)
    real(8),parameter::CFL=1000.0               !CFL number,implicit, 10E4~10E6
    real(8),parameter::eps=1.0E-10              !the converg_SAence precision  of the averag_SAe density
    
    !-----------------------------------------
    ! Roe scheme
    ! limiter
     real(8)::K = 5.0
    !-----------------------------
    !LU-SGS
    real(8)::omg = 1.5  !overrelaxation factor , 1< omg_SA <= 2
    !-------------------------------------------
    ! the scheme to calculate the properities on the edge, which is used in turbulent model.
    integer:: schemeType = 1                ! 1. one order upwind shceme
                                            ! 2. central difference scheme
     
    !---------------------------------------------------------------
    character(len=100)::caseRoute = "D:\Users\Hang\Documents\Visual Studio\Case"
    
    character(len=1),parameter:: gridNum = "1"             ! 1. triangle grid
                                                           ! 2. hybird grid                                                           
    character(len=1),parameter:: isNewGrid="y"             ! determine wheather new gird is used.
    
    !------------------------------------------------------------------------------
    character(len=1),parameter:: schemeNum = "3"            !  1. JST and SA
                                                            !  2. Roe and SA
                                                            !  3. JST and SST
                                                            !  4. Roe and SST
         
end module