subroutine artDissipation_SA( U,alf)                !artificial dissipation calculation
    implicit none
    real(8),intent(in)::U(:,:),alf(:)
    integer::i
    integer::ncl,ncr
    real(8)::dis                            !the artificial dissipation of every edge
    real(8)::nu,epsi2,epsi4                    !the pressure sensor
    real(8),allocatable::d2WT(:)                !the  median to calculate forth-order difference
    
    !write(*,*) "Art_dissipation"
    
    allocate( d2WT(ncells) )
    
    !initialize
    Dissi_SA = 0.0
    d2WT = 0.0
   
    !the first cycle to calculate the d2WT
    do i = 1,nedges
        ncl = iedge(3,i)                      
        ncr = iedge(4,i)
        if( ncr .GT. 0 )  then         !edge has no right cell is ignored
            d2WT(ncl) = d2WT(ncl) + 0.5*( W_SA(ncr) - W_SA(ncl) )
            d2WT(ncr) = d2WT(ncr) + 0.5*( W_SA(ncl) - W_SA(ncr) )
        end if
    end do
     
    !the second cycle to calculate the dissipation
    do i = 1,nedges
        ncl = iedge(3,i)
        ncr = iedge(4,i)
        
        if( ncr .GT. 0 )  then    !ncr less or equal to zero,the dissipation of this edge is set to zero
           
            nu = abs( U(5,ncr)-U(5,ncl) ) / abs( U(5,ncr)+U(5,ncl) ) 
            epsi2 = k2*nu
            epsi4 = max( 0.0 , k4-epsi2 )
            dis = alf(i)*epsi2*( W_SA(ncr) - W_SA(ncl) ) - alf(i)*epsi4*( d2WT(ncr) - d2WT(ncl) )
        
            !the left cell plus the dissipation and the right is contrary
            Dissi_SA(ncl) = Dissi_SA(ncl) + dis
            Dissi_SA(ncr) = Dissi_SA(ncr) - dis 
           
        end if
        
    end do
    
end subroutine
