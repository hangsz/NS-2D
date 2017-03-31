subroutine conFlux_JST      
    ! purpose:  use Jamson central difference scheme to calculate the convective flux.
    
    implicit none 
    integer::i
    integer::ncl,ncr                                  ! the number of the left and the right cell
    real(8)::rou0,c0                                  ! the reference value to calculate boundaries'S_SA flow information
    real(8)::c_av              
    real(8)::W_av(5)                                  ! conservative variables on any edge
    real(8)::ut,vt,Vnt,ux,uy,Vn                       ! define as U(2,:)*dy+U(3,:)*dx
    real(8)::Ma                                       ! mach number on any edge
    !-------------------------------------
    
    !write(*,*)  "Convective_flux"
    
    Fc = 0.0
    lamda_c = 0.0
    lamda_v = 0.0
    
    do i = 1,nedges
        ncl=iedge(3,i)
        ncr=iedge(4,i)
        
        ut =  0.5 * ( U_Rot(1,iedge(1,i)) + U_Rot(1,iedge(2,i)) ) 
        vt =  0.5 * ( U_Rot(2,iedge(1,i)) + U_Rot(2,iedge(2,i)) )
        Vnt = ut* vector(1,i)   + vt * vector(2,i) 
       
        ux =  U_av(2,i) - ut
        uy =  U_av(3,i) - vt
        Vn =  ux* vector(1,i) + uy * vector(2,i) 
        !Vn= U_av(2,i)*vector(1,i) + U_av(3,i)*vector(2,i)
        
        W_av(1)=U_av(1,i)
        W_av(2)=U_av(1,i)*U_av(2,i)
        W_av(3)=U_av(1,i)*U_av(3,i)
        W_av(5)=U_av(5,i)/(gamma-1) + U_av(1,i)*( U_av(2,i)**2+U_av(3,i)**2 )/2.0    !rouE 
        W_av(5)=W_av(5)+U_av(5,i)   !rouH  
        !-------------------------------------
        !propagation speed on each edge
        c_av=sqrt( gamma*R*U_av(6,i) )
        alf(i)= abs(Vn) + c_av * ds(i) 
        alf_v(i)= max( 4.0/3/U_av(1,i),gamma/U_av(1,i) ) *( muL_av(i)/PrL+muT_av(i)/PrT )
        
        !-------------------------------------
        
         !the left cell plus the convective flux
        Fc(1,ncl)= Fc(1,ncl) + Vn*W_av(1)
        Fc(2,ncl)= Fc(2,ncl) + Vn*W_av(2) + U_av(5,i)*vector(1,i)
        Fc(3,ncl)= Fc(3,ncl) + Vn*W_av(3) + U_av(5,i)*vector(2,i)
        Fc(5,ncl)= Fc(5,ncl) + Vn*W_av(5)
        Fc(5,ncl)= Fc(5,ncl) + Vnt * U_av(5,i)   
            
        lamda_c(ncl) = lamda_c(ncl) + alf(i)
        lamda_v(ncl) = lamda_v(ncl) + alf_v(i)*ds(i)**2/vol(ncl)
            
        if ( ncr .GT. 0) then
          
            !the right cell subtract the convective flux
            Fc(1,ncr)= Fc(1,ncr) - Vn*W_av(1)
            Fc(2,ncr)= Fc(2,ncr) -( Vn*W_av(2) + U_av(5,i)*vector(1,i) )
            Fc(3,ncr)= Fc(3,ncr) -( Vn*W_av(3) + U_av(5,i)*vector(2,i) )
            Fc(5,ncr)= Fc(5,ncr) - Vn*W_av(5)
            Fc(5,ncr)= Fc(5,ncr) - Vnt * U_av(5,i)    
            
            lamda_c(ncr) = lamda_c(ncr) + alf(i)
            lamda_v(ncr) = lamda_v(ncr) + alf_v(i)*ds(i)**2/vol(ncr)      
        end if        
    end do
    
    dt = CFL*vol/lamda_c 
   
end subroutine

    