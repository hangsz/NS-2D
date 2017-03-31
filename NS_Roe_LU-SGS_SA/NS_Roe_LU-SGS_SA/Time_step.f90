subroutine Time_step
    implicit none
    integer::i
    integer::ncl,ncr
  
    real(8)::Vn,c_av,mu_av!,alf_c,alf_v
    
    
    lamda_c = 0.0
    lamda_v = 0.0
    
    
    do i=1,nedges
        
        ncl=iedge(3,i)
        ncr=iedge(4,i)

        Vn = U_av(2,i)*vector(1,i) + U_av(3,i)*vector(2,i)
       
        c_av=sqrt( U_av(5,i)*gamma/U_av(1,i) )
        
        alf(i)= abs(Vn) + c_av * ds(i)
        alf_v(i)= max( 4.0/3/U_av(1,i),gamma/U_av(1,i) ) *( muL_av(i)/PrL+muT_av(i)/PrT )*ds(i)**2
        
        if(ncr .GT. 0) then
            lamda_c(ncl)  = lamda_c(ncl) + alf(i)
            lamda_c(ncr)  = lamda_c(ncr) + alf(i)
            
            lamda_v(ncl)  = lamda_v(ncl) + alf_v(i)/vol(ncl)
            lamda_v(ncr)  = lamda_v(ncr) + alf_v(i)/vol(ncr)
        else
            lamda_c(ncl)  = lamda_c(ncl) + alf(i)
            
            lamda_v(ncl)  = lamda_v(ncl) + alf_v(i)/vol(ncl)
        end if
        
    end do
    
    dt = CFL*vol/(lamda_c + C_time*lamda_v )   
    
end subroutine
