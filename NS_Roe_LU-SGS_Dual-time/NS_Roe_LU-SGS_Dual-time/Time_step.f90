subroutine Time_step
    implicit none
    integer::i
    integer::ncl,ncr
    real(8)::Vn,c_av,mu_av,alf_c,alf_v
    
    lamda_c = 0.0
    lamda_v = 0.0
     
    do i=1,nedges
        
        ncl=iedge(3,i)
        ncr=iedge(4,i)
        
        Vn = U_av(2,i)*vector(1,i) + U_av(3,i)*vector(2,i)
        c_av=sqrt( gamma*U_av(5,i)/U_av(1,i) )
        
       
        alf_c= abs(Vn) + c_av * ds(i)
        
        mu_av = 1.45*U_av(6,i)**(3.0/2)/(U_av(6,i)+110.0) *1.0E-6
        alf_v= max( 4.0/3/U_av(1,i),gamma/U_av(1,i) ) *( mu_av/PrL)*ds(i)**2
        
        if(ncr .GT. 0) then
            lamda_c(ncl)  = lamda_c(ncl) + alf_c
            lamda_c(ncr)  = lamda_c(ncr) + alf_c
            
            lamda_v(ncl)  = lamda_v(ncl) + alf_v/vol(ncl)
            lamda_v(ncr)  = lamda_v(ncr) + alf_v/vol(ncr)
        else
            lamda_c(ncl)  = lamda_c(ncl) + alf_c
            
            lamda_v(ncl)  = lamda_v(ncl) + alf_v/vol(ncl)
        end if
        
    end do
    
    dt = CFL*vol/(lamda_c + C_time*lamda_v )   
    
end subroutine
