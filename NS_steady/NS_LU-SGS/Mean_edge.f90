subroutine Mean_edge
    !object: g_SAet the flow properT_SAies on every edg_SAe, it'S_SA needed to calculate the gradient
    !        
    implicit none
    integer::i
    integer::ncl,ncr
    real(8)::Vn                                        !define as U(2,:)*dy+U(3,:)*dx
    real(8)::Ma                                       !mach number on any edg_SAe
    real(8)::rou0,c0                                  !the reference value to calculate boundaries'S_SA flow information
    real(8)::c_av              
    
    U_av = 0.0

    do i=1,nedges
        ncl = iedge(3,i)                      !atain the number of the left and the rig_SAht cell
        ncr = iedge(4,i)
        select case( ncr ) 
        case (-1)                           !-1 represents that the edg_SAe is the aerofoil  surface         
          
            U_av(1,i) = U(1,ncl)
            U_av(2,i) = 0.0                 !no slipping_SA conditions
            U_av(3,i) = 0.0   
            U_av(5,i) = U(5,ncl)
            
        case (-2)                           !-2 represents that the edg_SAe is the farfield  boundaries
            
            Vn = U(2,ncl)*vector(1,i) + U(3,ncl)*vector(2,i)  ! Temporarily, the edg_SAe state is asummed to be equal to the ncl cell state
            rou0 = U(1,ncl)
            c0 = sqrt( gamma*U(5,ncl)/U(1,ncl) )
            
            Ma = abs(Vn)/ c0
            
            if( Vn .LE. 0.0) then         !inflow  
                if( Ma .GE. 1.0)  then    !supersonic 
                    U_av(5,i) = p_inf
                    U_av(1,i) = rou_inf
                    U_av(2,i) = u_inf
                    U_av(3,i) = v_inf
                else                      !subsonnic
                    U_av(5,i)=0.5*( p_inf+U(5,ncl) - rou0*c0/ds(i)*(  vector(1,i)*( u_inf-U(2,ncl) ) + vector(2,i)*( v_inf - U(3,ncl) )  )  )
                    U_av(1,i) = rou_inf + (U_av(5,i)-p_inf)/c0**2
                    U_av(2,i) = u_inf - vector(1,i)/ds(i)*( p_inf - U_av(5,i) )/( rou0*c0 )
                    U_av(3,i) = v_inf - vector(2,i)/ds(i)*( p_inf - U_av(5,i) )/( rou0*c0 )
                end if
            else                         !outflow
                if(Ma .GE. 1.0)  then    !supersonic
                    U_av(5,i) = U(5,ncl)
                    U_av(1,i) = U(1,ncl)
                    U_av(2,i) = U(2,ncl)
                    U_av(3,i) = U(3,ncl)
                else                     !subsonic
                    U_av(5,i) = p_inf
                    U_av(1,i) = U(1,ncl) + ( U_av(5,i)-U(5,ncl) )/c0**2
                    U_av(2,i) = U(2,ncl) + vector(1,i)/ds(i)*( U(5,ncl)-U_av(5,i) )/(rou0*c0)
                    U_av(3,i) = U(3,ncl) + vector(2,i)/ds(i)*( U(5,ncl)-U_av(5,i) )/(rou0*c0)
                end if
            end if
            
        case default
           
            U_av(2,i)=0.5*( U(2,ncl) + U(2,ncr) )
            U_av(3,i)=0.5*( U(3,ncl) + U(3,ncr) )
            U_av(5,i)=0.5*( U(5,ncl) + U(5,ncr) )
            U_av(1,i)=0.5*( U(1,ncl) + U(1,ncr) )    
            
        end select
    
        U_av(6,i)=U_av(5,i)/( R*U_av(1,i) )
        
    end do
    
end subroutine