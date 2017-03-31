subroUTine Con_flux_SST( U_av)       !convective flux
    implicit none 
    real(8),intent(in)::U_av(:,:)     
    
    integer::i
    integer::ncl,ncr                                       !the number of the left and the right cell
                                      
    real(8)::WT_av(2)
    real(8)::Vn                                     
                              
    !set Fc_SST to zero
   
    Fc_SST = 0.0
    do i = 1,nedges
        ncl=iedge(3,i)
        ncr=iedge(4,i)
        
        Vn= U_av(2,i)*vector(1,i) + U_av(3,i)*vector(2,i)
       
        WT_av = U_av(1,i)*U_av_SST(:,i)
       
        if ( ncr .GT. 0) then
            !the left cell plus the convective flux
            
            Fc_SST(:,ncl)  = Fc_SST(:,ncl) + WT_av * Vn
            Fc_SST(:,ncr)  = Fc_SST(:,ncr) - WT_av * Vn
        else 
            Fc_SST(:,ncl)  = Fc_SST(:,ncl) + WT_av * Vn
        end if

    end do
    
end subroUTine
