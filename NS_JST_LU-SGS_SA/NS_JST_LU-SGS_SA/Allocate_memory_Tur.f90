subroutine Allocate_memory_Tur

    implicit none
    
    allocate( nuT(ncells) )
    
    allocate( WT(ncells) )
    
    allocate( nuT_av(nedges) )
    allocate( muT(ncells) ) 
    allocate( muT_av(nedges) )
    
    allocate( FcT(ncells) ) 
    
    allocate( Grad_nuT(2,ncells) )
    allocate( Grad_WT(2,ncells) )
    
    allocate( FvT(ncells) ) 
    allocate( Dissi_T(ncells) )
    allocate( QT(ncells) )
    allocate( dQTdWT(ncells) )
    
    allocate( RsiT(ncells) )
    
    
    allocate( S(ncells) )
    allocate( chi(ncells) )
    allocate( fv1(ncells) )
    allocate( fv2(ncells) )
    allocate( fv3(ncells) )
    
    allocate( fw(ncells) )
    allocate( g(ncells) )
    allocate( rT(ncells) )
    
end subroutine
