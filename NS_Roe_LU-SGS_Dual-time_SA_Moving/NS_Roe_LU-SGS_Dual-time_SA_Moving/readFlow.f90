subroutine readFlow  
    !----------------------------------------------------
        ! purpose: read the stady flow properties
    !----------------------------------------------------
        
    implicit none
    integer::i
    integer::idRead
    character( len =30)::filename
    logical::alive
    
    !----------------------------------------------------
    
    filename = 'steady/steady.dat'
    inquire(file =filename,exist =alive)
    
    if ( .NOT. alive) then
        write(*,*)  filename,'not exist.'
        stop   
    end if
    
    open( newunit=idRead,file=filename)
    
    read(idRead,*) 
    read(idRead,*) 
    read(idRead,*) 
    read(idRead,*) 
    read(idRead,*) 
     
    do i=1,nnodes
        read(idRead,*) !xy(1,i)
    end do
    do i=1,nnodes
        read(idRead,*) !xy(2,i)
    end do
    
    do i=1,ncells
        read(idRead,*) U(5,i)
    end do
    do i=1,ncells
        read(idRead,*) U(1,i)
    end do
    do i=1,ncells
        read(idRead,*) U(2,i)
    end do
    do i=1,ncells
        read(idRead,*) U(3,i)
    end do
    do i=1,ncells
        read(idRead,*) nuT(i)
    end do
    do i=1,ncells
        read(idRead,*) !icell(1,i),icell(2,i),icell(3,i)
    end do
    
    close(idRead)
   
end subroutine