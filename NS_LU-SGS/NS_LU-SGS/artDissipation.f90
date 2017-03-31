subroutine artDissipation                !artificial dissipation calculation
    implicit none
    integer::i
    integer::ncl,ncr
    real(8)::dis(5)                            !the arTificial dissipation of every edge
    real(8)::nu,epsi2,epsi4                    !the pressure sensor
    real(8),allocatable::d2W(:,:)                !the  median to calculate forth-order difference
    
    !write(*,*) "Art_dissipation"
    allocate(d2W(5,ncells))
    !initialize
    Dissi = 0.0
    d2W = 0.0 
   
    !the first cycle to calculate the d2W
    do i = 1,nedges
        ncl = iedge(3,i)                      
        ncr = iedge(4,i)
        if( ncr .GT. 0 )  then         !edge has no right cell is ignored
            !d2W(:,ncl) = d2W(:,ncl) + 0.5*( W(:,ncr) - W(:,ncl) )
            !d2W(:,ncr) = d2W(:,ncr) + 0.5*( W(:,ncl) - W(:,ncr) )
            d2W(:,ncl) = d2W(:,ncl) + ( W(:,ncr) - W(:,ncl) )
            d2W(:,ncr) = d2W(:,ncr) + ( W(:,ncl) - W(:,ncr) )
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
            dis = alf(i)*epsi2*( W(:,ncr) - W(:,ncl) ) - alf(i)*epsi4*( d2W(:,ncr) - d2W(:,ncl) )
            !the sum of every edge'S_SA spectral radius
               
            
            !the left cell plus the dissipation and the right is contrary
            Dissi(:,ncl) = Dissi(:,ncl) + dis
            Dissi(:,ncr) = Dissi(:,ncr) - dis 
                
        end if
        
    end do
    
end subroutine
