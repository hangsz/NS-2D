subroutine meanEdge  
    !------------------------------------------
        ! purpose: get the flow properties on each edge, it's necessary to calculate the gradient
    !------------------------------------------        
    implicit none
    integer::i
    integer::ncl,ncr
    real(8)::ux,uy,Vn                   ! the flow speed cross each edge
    real(8)::Ma                         ! mach number on any edge 
    real(8)::rou0,c0                    ! the reference value to calculate boundaries's flow information
    real(8)::c_av                       ! the sound speed of each edge
    
    U_av = 0.0                                    
    
    do i=1,nedges
        ncl = iedge(3,i)                    !get the number of the left and the right cell
        ncr = iedge(4,i)
        select case( ncr ) 
        case (-1)                           ! -1 represents that the edge is the wall, using wall boundary conditions
          
            U_av(1,i) = U(1,ncl)            ! drou/dy = 0, normal gradient of density is zero
            U_av(2,i) = 0.5 * ( U_Rot(1,iedge(1,i)) + U_Rot(1,iedge(2,i)) )                 ! no slipping condition
            U_av(3,i) = 0.5 * ( U_Rot(2,iedge(1,i)) + U_Rot(2,iedge(2,i)) )   
            U_av(5,i) = U(5,ncl)            ! dP/dt = 0, normal gradient of pressure is zero
            
        case (-2)                           ! -2 represents that the edge is the farfield  boundaries, using farfield boundary conditions
            
            ! the real speed of the flow crossing the edge has to subtract the grid speed from the folw speed
            ux =  U(2,ncl) - 0.5 * ( U_Rot(1,iedge(1,i)) + U_Rot(1,iedge(2,i)) )
            uy =  U(3,ncl) - 0.5 * ( U_Rot(2,iedge(1,i)) + U_Rot(2,iedge(2,i)) )
            Vn = ux* vector(1,i)   + uy * vector(2,i)
           ! Vn = U(2,ncl)*vector(1,i) + U(3,ncl)*vector(2,i)  
           
           ! Temporarily, the edge state is asummed to be equal to the ncl cell state
            rou0 = U(1,ncl)
            c0 = sqrt( gamma*U(5,ncl)/U(1,ncl) )
            
            Ma = abs(Vn) / c0
            
            !Ma = sqrt(U(2,ncl)**2 + U(3,ncl)**2) / c0
            
            if( Vn .LE. 0.0) then         ! inflow  
                if( Ma .GE. 1.0)  then    
                    ! supersonic 
                    ! convecition dominated, flow properities on the inflow boundaries are determined by the inflow 
                    U_av(5,i) = p_inf    
                    U_av(1,i) = rou_inf
                    U_av(2,i) = u_inf
                    U_av(3,i) = v_inf
                else                     
                    ! subsonnic
                    ! convection does't dominated, flow is infuenced by both the upwind flow and the downwind flow
                    U_av(5,i)=0.5*( p_inf+U(5,ncl) - rou0*c0/ds(i)*(  vector(1,i)*( u_inf-U(2,ncl) ) + vector(2,i)*( v_inf - U(3,ncl) )  )  )
                    U_av(1,i) = rou_inf + (U_av(5,i)-p_inf)/c0**2
                    U_av(2,i) = u_inf - vector(1,i)/ds(i)*( p_inf - U_av(5,i) )/( rou0*c0 )
                    U_av(3,i) = v_inf - vector(2,i)/ds(i)*( p_inf - U_av(5,i) )/( rou0*c0 )
                    
                end if
                
            else                          ! outflow
                if(Ma .GE. 1.0)  then    
                    ! supersonic
                    ! convection dominated, flow propertities on the outer boundaries are determined by the outflow
                    U_av(5,i) = U(5,ncl)  
                    U_av(1,i) = U(1,ncl)
                    U_av(2,i) = U(2,ncl)
                    U_av(3,i) = U(3,ncl)
                else                    
                    ! subsonic
                    ! convection does't dominated, flow is infuenced by both the upwind flow and the downwind flow
                    U_av(5,i) = p_inf
                    U_av(1,i) = U(1,ncl) + ( U_av(5,i)-U(5,ncl) )/c0**2
                    U_av(2,i) = U(2,ncl) + vector(1,i)/ds(i)*( U(5,ncl)-U_av(5,i) )/(rou0*c0)
                    U_av(3,i) = U(3,ncl) + vector(2,i)/ds(i)*( U(5,ncl)-U_av(5,i) )/(rou0*c0)
                    
                end if
                
            end if
            
        case default       ! folw properities are averaged on the edge in the internal domain        
           
            U_av(2,i)=0.5*( U(2,ncl) + U(2,ncr) )
            U_av(3,i)=0.5*( U(3,ncl) + U(3,ncr) )
            U_av(5,i)=0.5*( U(5,ncl) + U(5,ncr) )
            U_av(1,i)=0.5*( U(1,ncl) + U(1,ncr) )    
            
        end select
       
        ! temperature are calculated by the ideal gas law
        U_av(6,i)=U_av(5,i)/( R*U_av(1,i) )
        
    end do
    
end subroutine