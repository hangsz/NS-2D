subroUTine conFlux_SST( U_av,U_Rot)       !convective flux
    implicit none 
    real(8),intent(in)::U_av(:,:),U_Rot(:,:)     
    
    integer::i
    integer::ncl,ncr                                       !the number of the left and the right cell
                                      
    real(8)::WT_av(2)
    real(8)::ux,uy,Vn     
    
                              
    lamda_SST = 0.0
    
    Fc_SST = 0.0
    do i = 1,nedges
        ncl=iedge(3,i)
        ncr=iedge(4,i)
        
        ux =  U_av(2,i) - 0.5 * ( U_Rot(1,iedge(1,i)) + U_Rot(1,iedge(2,i)) )
        uy =  U_av(3,i) - 0.5 * ( U_Rot(2,iedge(1,i)) + U_Rot(2,iedge(2,i)) )
        Vn = ux* vector(1,i)   + uy * vector(2,i)
        !Vn= U_av(2,i)*vector(1,i) + U_av(3,i)*vector(2,i)
       
        WT_av = U_av(1,i)*U_av_SST(:,i)
        
           !-------------------------------------
        !propagation speed on each edge
        alf_SST(i)= abs(Vn)
        
        ! alf_v_T(i) = 0.0 !alf_v(i)
        ! alf_v_T(i)= max( 4.0/3/U_av(1,i),gamma/U_av(1,i) ) *( muL_av(i)/PrL+muT_av(i)/PrT )*ds(i)**2   !???????????????????????????????
        
        Fc_SST(:,ncl)  = Fc_SST(:,ncl) + WT_av * Vn
            
        lamda_SST(ncl) = lamda_SST(ncl) + alf_SST(i)
        if ( ncr .GT. 0) then
            !the left cell plus the convective flux
            
            Fc_SST(:,ncr)  = Fc_SST(:,ncr) - WT_av * Vn
            lamda_SST(ncr) = lamda_SST(ncr) + alf_SST(i)
             
            !!lamda_v_T(ncl) = lamda_v_T(ncl) + alf_v_T(i)/vol(ncl)
            !!lamda_v_T(ncr) = lamda_v_T(ncr) + alf_v_T(i)/vol(ncr)
            !!
        else 
            
            !lamda_v_T(ncl) = lamda_v_T(ncl) + alf_v_T(i)/vol(ncl)

        end if

    end do
    
end subroUTine
