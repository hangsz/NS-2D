subroutine meanEdge  
    !------------------------------------------
        ! purpose: get the flow properties on each edge, it's necessary to calculate the gradient 
        !          using Rimenna invarient variables
    !------------------------------------------        
    implicit none
    integer::i
    integer::ncl,ncr
    real(8)::nx,ny
    real(8)::ux,uy,Vn                   ! the flow speed cross each edge
    real(8)::Ma                         ! mach number on any edge 
    real(8)::R_p,R_n,Vnor,Vtau,s,c! the sound speed of each edge
    
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
            
            nx = vector(1,i)/ds(i)
            ny = vector(2,i)/ds(i)
            ! the real speed of the flow crossing the edge has to subtract the grid speed from the folw speed
            ux =  U(2,ncl) - 0.5 * ( U_Rot(1,iedge(1,i)) + U_Rot(1,iedge(2,i)) )
            uy =  U(3,ncl) - 0.5 * ( U_Rot(2,iedge(1,i)) + U_Rot(2,iedge(2,i)) )
            Vn = ux * nx + uy * ny
            
            ! Temporarily, the edge state is asummed to be equal to the ncl cell state
            Ma = abs(Vn) / sqrt( gamma*U(5,ncl)/U(1,ncl) )
            
            if( Vn .LE. 0.0) then         ! inflow  
                if( Ma .GE. 1.0)  then    
                    ! supersonic 
                    ! convecition dominated, flow properities on the inflow boundaries are determined by the inflow 
                    
                    U_av(2,i) = u_inf
                    U_av(3,i) = v_inf
                    U_av(1,i) = rou_inf
                    U_av(5,i) = p_inf    
                else                     
                    ! subsonnic
                    ! convection does't dominated, flow is infuenced by both the upwind flow and the downwind flow
                    R_p =  U(2,ncl)*nx + U(3,ncl)*ny + 2*sqrt( gamma*U(5,ncl)/U(1,ncl) )/( gamma-1)                                           ! R+
                    R_n =  u_inf*nx + v_inf*ny - 2*a_inf/(gamma-1)
                   
                    Vnor = 0.5*( R_p + R_n)
                    c = (gamma-1)/4*( R_p - R_n )
                    Vtau = u_inf*ny - v_inf*nx   !  edge node1 => node2
                    s = p_inf/rou_inf**gamma
                    
                    U_av(2,i) = Vnor * nx + Vtau * ny
                    U_av(3,i) = Vnor * ny - Vtau * nx
                    U_av(1,i) = ( c**2/( gamma*s) )**( 1.0/(gamma-1) )
                    U_av(5,i) = s*U_av(1,i)**gamma
                       
                end if
                
            else    
                ! outflow
                if(Ma .GE. 1.0)  then    
                    ! supersonic
                    ! convection dominated, flow propertities on the outer boundaries are determined by the outflow
                    
                    U_av(2,i) = U(2,ncl)
                    U_av(3,i) = U(3,ncl)
                    U_av(1,i) = U(1,ncl)
                    U_av(5,i) = U(5,ncl)  
                else  
                    
                    R_p =  U(2,ncl)*nx + U(3,ncl)*ny + 2*sqrt( gamma*U(5,ncl)/U(1,ncl) )/( gamma-1)                                           ! R+
                    R_n =  u_inf*nx + v_inf*ny - 2*a_inf/(gamma-1)
                    
                    Vnor = 0.5*( R_p + R_n)
                    c = (gamma-1)/4*( R_p - R_n )
                    Vtau = U(2,ncl)*ny - U(3,ncl)*nx   !  edge node1 => nodes2
                    s = U(5,ncl)/U(1,ncl)**gamma
                    
                    U_av(2,i) = Vnor*nx + Vtau *ny
                    U_av(3,i) = Vnor*ny - Vtau *nx
                    U_av(1,i) = ( c**2/( gamma*s) )**( 1.0/(gamma-1) )
                    U_av(5,i) = s*U_av(1,i)**gamma
                        
                end if
                
            end if
            
        case default       ! folw properities are averaged on the edge in the internal domain        
           
            U_av(2,i)=0.5*( U(2,ncl) + U(2,ncr) )
            U_av(3,i)=0.5*( U(3,ncl) + U(3,ncr) )
            U_av(1,i)=0.5*( U(1,ncl) + U(1,ncr) )    
            U_av(5,i)=0.5*( U(5,ncl) + U(5,ncr) )
            
        end select
       
        ! temperature are calculated by the ideal gas law
        U_av(6,i)=U_av(5,i)/( R*U_av(1,i) )
        
    end do
    
    muL =   1.45*( U(5,:)/R/U(1,:) )**(3.0/2)/ ( U(5,:)/R/U(1,:) + 110.0 ) * 1.0E-6
    muL_av = 1.45* U_av(6,:)**(3.0/2)/( U_av(6,:)+110.0 ) *1.0E-6
       
end subroutine