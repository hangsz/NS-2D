module  controlPara    !Control parameters module
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
    !Finite volume method 
    
    real(8),parameter::k2=0.6,k4=1.0/64         !arT_SAifical dissipation,k2(1/2,1),k4(1/256,1/32)
    !real(8),parameter::CFL=2.0                 !CFL number,explicit, less than 2*sqrt(2)
    real(8),parameter::CFL=10000.0               !CFL number,implicit, 10E4~10E6
    real(8),parameter::eps=1.0E-4              !the converg_SAence precision  of the averag_SAe density
    
   
    !-------------------------------------------
    !LU-SGS,is used to calculate the D matrix
    real(8),parameter::omg = 1.5  !overrelaxation factor , 1< omg <= 2
    
    !-------------------------------------------
    ! the scheme to calculate the properities on the edge, which is used in turbulent model.
    integer,parameter:: schemeType = 2               
    !                   1. one order upwind shceme                             
    !                   2. central difference scheme                                    
    !------------------------------------------------------
                                            
    !dual time
    integer,parameter::iter_inner = 500           ! the dual time scheme iterative times 
    real(8),parameter::eps_inner = 1.0E-3
    !-----------------------------------------------------------------
    !                 static  or dynamic                             !
    integer,parameter::moveNum = 1                                   !
    !                    0.  static                                  !
    !                    1.  rotate around the pressure center       !
    !-----------------------------------------------------------------
    integer,parameter::iterMax=10000           ! static 
    
    integer,parameter::stepMax=50              ! steps every 1/4 period
    integer,parameter::steps = 25              ! if reach the steps, output   
    
    integer,parameter::phase= 24                       ! the phase, three period( 2 phase every 1/4 period )
    !real(8),parameter::dt_r = 0.00095                  ! real time step 
    real(8),parameter::rot_cen(2)=[0.273,0.0]          !rotational center,(m,m)
    real(8),parameter::att_ampl = 2.41                 !rotational amplitude,бу
    real(8),parameter::kr = 0.0808                     !reduced frequency
    
    !-----------------------------------------------------------------
    !              onvective flux scheme                             !
    integer,parameter::conFluxNum = 1                                !      
    !                   1.  JST scheme                               !       
    !                   2.  Roe scheme                               !
    !                                                                !
    ! Roe scheme                                                     !
    ! limiter                                                        !
    real(8),parameter::K = 5.0                                       !
    !-----------------------------------------------------------------
    !-----------------------------------------------------------------
    !                 turbulent model  select                        !
    integer,parameter:: turModelNum = 0                              !            
    !                       0.  doesn't use turbulent model            !
    !                       1.  SA                                   !              
    !                       2.  SST                                  !                                   
    !-----------------------------------------------------------------
    !                     set the case route and grid number
    character(len=100)::caseRoute = "D:\Users\Hang\Documents\Visual Studio\Case"   ! case route
    !
    character(len=3),parameter:: gridNum = "2"             ! 1. triangle grid
    !                                                      ! 2. hybird grid                                                           
    character(len=1),parameter:: isNewGrid="y"             !    determine wheather new gird is used.
    !
    !------------------------------------------------------------------------------
    !                  determine which file as the word directory                 !
    character(len=3),parameter:: fileNum = "1"                                    !
    !                       1. JST                                                !
    !                       2. Roe                                                !
    !                       3. JST and SA                                         !
    !                       4. Roe and SA                                         !
    !                       5. JST and SST                                        !
    !                       6. Roe and SST                                        !
    !------------------------------------------------------------------------------     
    end module
    
    
    !-------------------------------------------------------------------------------
    !                   test  specification             
    ! grid number              1                     2                       3
    ! JST                    succeed               fail
    ! Roe                    succeed               fail
    ! JST SA                 succeed               succeed                               
    ! Roe SA                 succeed               succeed                                   
    ! JST kw-SST             succeed               fail                                      
    ! Roe kw-SST              -                       -                             
    !-------------------------------------------------------------------------------
    
    
    !-----------------------------------------------------------------------
    ! version 1.0                10-27-14               modify the fifth convetive flux  + Vt*p