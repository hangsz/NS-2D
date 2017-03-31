subroutine Read_flow
    implicit none
    integer::i
    integer::fileid = 7
    
    open(fileid,file="steady.dat")
    
    read(fileid,*) 
    read(fileid,*) 
    read(fileid,*) 
    read(fileid,*) 
    read(fileid,*) 
     
    do i=1,nnodes
        read(fileid,*) !xy(1,i)
    end do
    do i=1,nnodes
        read(fileid,*) !xy(2,i)
    end do
    
    do i=1,ncells
        read(fileid,*) U(5,i)
    end do
    do i=1,ncells
        read(fileid,*) U(1,i)
    end do
    do i=1,ncells
        read(fileid,*) U(2,i)
    end do
    do i=1,ncells
        read(fileid,*) U(3,i)
    end do
    do i=1,ncells
        read(fileid,*) UT(1,i)
    end do
    do i=1,ncells
        read(fileid,*) UT(2,i)
    end do
    do i=1,ncells
        read(fileid,*) !icell(1,i),icell(2,i),icell(3,i)
    end do
    
    close(fileid)
   
end subroutine