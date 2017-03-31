subroutine Output(filename)
    implicit none 
    integer::i
    integer,parameter::fileid=7
    character(len =30)::filename
    !write(*,*)   "Output"
    
    
    !********************************************
    open(fileid,file="flow/"//trim(filename)//".dat")
    
    write(fileid,*) 'TITLE="euler solver on unstructured grid" '
    write(fileid,*) 'VARIABLES= "X" "Y"  "P" "Rou" "U"  "V" "nuT"' 
    write(fileid,*) "ZONE NODES=",nnodes," ELEMENTS=",ncells," ZONETYPE=FETRIANGLE"
    write(fileid,*) "DATAPACKING=BLOCK"
    write(fileid,*) "VARLOCATION=([3-7]=CELLCENTERED)"
     
    do i=1,nnodes
        write(fileid,*) xy(1,i)
    end do
    do i=1,nnodes
        write(fileid,*) xy(2,i)
    end do
    
    do i=1,ncells
        write(fileid,*) U(5,i)
    end do
    do i=1,ncells
        write(fileid,*) U(1,i)
    end do
    do i=1,ncells
        write(fileid,*) U(2,i)
    end do
    do i=1,ncells
        write(fileid,*) U(3,i)
    end do
    do i=1,ncells
        write(fileid,*) nuT(i)
    end do
    do i=1,ncells
        write(fileid,*) icell(1,i),icell(2,i),icell(3,i)
    end do
    
    close(fileid)

    
    

end subroutine 
