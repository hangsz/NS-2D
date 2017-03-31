subroutine Con_flux        !convective flux
    implicit none 
    integer::i
    integer::ncl,ncr                                  !the number of the left and the right cell
    real(8)::rou0,c0                                  !the reference value to calculate boundaries's flow information
  
    real(8)::c_av ,mu_av             
    real(8)::W_av(5)                                  !conservative variables on any edge
    real(8)::Vn                                       !define as U(2,:)*dy+U(3,:)*dx
    real(8)::Ma                                       !mach number on any edge
   
    !write(*,*)  "Convective_flux"
    
    !set Fc to zero
   
    Fc = 0.0
    do i = 1,nedges
        ncl=iedge(3,i)
        ncr=iedge(4,i)
       
      
        Vn= U_av(2,i)*vector(1,i) + U_av(3,i)*vector(2,i)
        W_av(1)=U_av(1,i)
        W_av(2)=U_av(1,i)*U_av(2,i)
        W_av(3)=U_av(1,i)*U_av(3,i)
        W_av(5)=U_av(5,i)/(gamma-1) + U_av(1,i)*( U_av(2,i)**2+U_av(3,i)**2 )/2.0    !rouE 
        W_av(5)=W_av(5)+U_av(5,i)   !rouH  
        
        if ( ncr .GT. 0) then
            !the left cell plus the convective flux
            Fc(1,ncl)= Fc(1,ncl) + Vn*W_av(1)
            Fc(2,ncl)= Fc(2,ncl) + Vn*W_av(2) + U_av(5,i)*vector(1,i)
            Fc(3,ncl)= Fc(3,ncl) + Vn*W_av(3) + U_av(5,i)*vector(2,i)
            Fc(5,ncl)= Fc(5,ncl) + Vn*W_av(5)
            !the right cell subtract the convective flux
            Fc(1,ncr)= Fc(1,ncr) - Vn*W_av(1)
            Fc(2,ncr)= Fc(2,ncr) -( Vn*W_av(2) + U_av(5,i)*vector(1,i) )
            Fc(3,ncr)= Fc(3,ncr) -( Vn*W_av(3) + U_av(5,i)*vector(2,i) )
            Fc(5,ncr)= Fc(5,ncr) - Vn*W_av(5)    
         
        else
            !the left cell plus the convective flux
            Fc(1,ncl)= Fc(1,ncl) + Vn*W_av(1)
            Fc(2,ncl)= Fc(2,ncl) + Vn*W_av(2) + U_av(5,i)*vector(1,i)
            Fc(3,ncl)= Fc(3,ncl) + Vn*W_av(3) + U_av(5,i)*vector(2,i)
            Fc(5,ncl)= Fc(5,ncl) + Vn*W_av(5)
        end if
        
        !time step
        !propagation speed on every edge
        c_av=sqrt( U_av(5,i)*gamma/U_av(1,i) )
        alf(i)= abs(Vn) + c_av * ds(i) 
        
        mu_av = 1.45*U_av(6,i)**(3.0/2)/(U_av(6,i)+110.0) *1.0E-6
        alf_v(i)= max( 4.0/3/U_av(1,i),gamma/U_av(1,i) ) *( mu_av/PrL )*ds(i)**2
    end do
    
end subroutine

    