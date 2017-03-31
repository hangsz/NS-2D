subroutine Allocate_memory_Tur

    implicit none
    
    allocate( muT(ncells) )
    allocate( muT_av(nedges) )
     
    allocate( UT(2,ncells) )
    allocate( UT_av(2,nedges) )
    
    allocate( WT(2,ncells) )
    
    allocate( WTn(2,ncells) )
    allocate( WTn1(2,ncells) )
    allocate( Q_WT(2,ncells) )
    
    
    allocate( FcT(2,ncells) ) 
    
    allocate( Grad_UT(2,2,ncells) )
    allocate( Grad_WT(2,2,ncells) )
    
    allocate( FvT(2,ncells) ) 
    allocate( Dissi_T(2,ncells) )
    allocate( QT(2,ncells) )
    allocate( dQTdWT(2,ncells) )
    
    allocate( RsiT(2,ncells) )
    
     allocate( f2(ncells) )   !LU-SGS
    
    allocate( sigma_K(ncells) )
    allocate( sigma_omg(ncells) )
    allocate( Cw(ncells) )
    allocate( beta(ncells) )
    
end subroutine
