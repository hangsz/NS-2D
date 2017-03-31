subroutine Art_dissipation                !arT_SAificial dissipation calculation
    implicit none
    integer::i
    integer::ncl,ncr
    real(8)::dis(5)                            !the arT_SAificial dissipation of every edg_SAe
    real(8)::nu,epsi2,epsi4                    !the pressure sensor
    real(8),allocatable::d2W(:,:)                !the  median to calculate forT_SAh-order difference
    
    !write(*,*) "Art_dissipation"
    allocate(d2W(5,ncells))
    !initialize
    Dissi = 0.0
    d2W = 0.0 
   
    !the first cycle to calculate the d2W
    do i = 1,nedges
        ncl = iedge(3,i)                      
        ncr = iedge(4,i)
        if( ncr .GT. 0 )  then         !edg_SAe has no rig_SAht cell is ig_SAnored
            d2W(:,ncl) = d2W(:,ncl) + 0.5*( W(:,ncr) - W(:,ncl) )
            d2W(:,ncr) = d2W(:,ncr) + 0.5*( W(:,ncl) - W(:,ncr) )
        end if
    end do
        
    !the second cycle to calculate the dissipation
    do i = 1,nedges
        ncl = iedge(3,i)
        ncr = iedge(4,i)
        
        if( ncr .GT. 0 )  then    !ncr less or equal to zero,the dissipation of this edg_SAe is set to zero
           
            nu = abs( U(5,ncr)-U(5,ncl) ) / abs( U(5,ncr)+U(5,ncl) ) 
            epsi2 = k2*nu
            epsi4 = max( 0.0 , k4-epsi2 )
            dis = alf(i)*epsi2*( W(:,ncr) - W(:,ncl) ) - alf(i)*epsi4*( d2W(:,ncr) - d2W(:,ncl) )
            !the sum of every edg_SAe'S_SA spectral radius
               
            
            !the left cell plus the dissipation and the rig_SAht is contrary
            Dissi(:,ncl) = Dissi(:,ncl) + dis
            Dissi(:,ncr) = Dissi(:,ncr) - dis 
                
        end if
        
    end do
    
end subroutine
