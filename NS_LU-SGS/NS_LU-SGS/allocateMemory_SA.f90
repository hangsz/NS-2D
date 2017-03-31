subroutine allocateMemory_SA

    write(*,*) "allocateMemory_SA"
    
    allocate( nu_SA(ncells) )
    
    allocate( W_SA(ncells) )
    !---------------
    ! dual time
    if( moveNum /= 0 ) then
        allocate( W_SA_n(ncells) )
        allocate( W_SA_n1(ncells) )
        allocate( QW_SA(ncells) )
    end if
    !--------------
    allocate( nu_av_SA(nedges) )
    
    allocate( Fc_SA(ncells) ) 
    
    allocate( Grad_nu_SA(2,ncells) )
    allocate( Grad_W_SA(2,ncells) )
    
    allocate( Fv_SA(ncells) ) 
    allocate( Dissi_SA(ncells) )
    allocate( Q_SA(ncells) )
    allocate( dQdW_SA(ncells) )
    
    allocate( Rsi_SA(ncells) )
    
    
    allocate( S_SA(ncells) )
    allocate( chi_SA(ncells) )
    allocate( fv1_SA(ncells) )
    allocate( fv2_SA(ncells) )
    allocate( fv3_SA(ncells) )
    
    allocate( fw_SA(ncells) )
    allocate( g_SA(ncells) )
    allocate( rT_SA(ncells) )
    
end subroutine
