subroutine allocateMemory_SST

    write(*,*) "allocateMemory_SST"
    
     
    allocate( U_SST(2,ncells) )
    allocate( U_av_SST(2,nedges) )
    
    allocate( W_SST(2,ncells) )
    
    !-----------
    ! dual time
    if( moveNum /= 0) then
        allocate( W_SST_n(2,ncells) )
        allocate( W_SST_n1(2,ncells) )
        allocate( QW_SST(2,ncells) )
    end if
    
    !-----------
    
    allocate( Fc_SST(2,ncells) ) 
    
    allocate( Grad_UT(2,2,ncells) )
    
    allocate( Fv_SST(2,ncells) ) 
    allocate( Dissi_SST(2,ncells) )
    allocate( Q_SST(2,ncells) )
    allocate( dQdW_SST(2,ncells) ) 
    
    allocate( Rsi_SST(2,ncells) )
    
     allocate( f2_SST(ncells) )       ! LU_SGS
    
    allocate( sigma_K_SST(ncells) )
    allocate( sigma_omg_SST(ncells) )
    allocate( Cw_SST(ncells) )
    allocate( beta_SST(ncells) )
    
end subroUTine
