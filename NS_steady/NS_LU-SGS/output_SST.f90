subroutine Output_SST(filename)
    implicit none 
    integer::i
    integer,parameter::fileid=7
    
    character(len=30)::filename
    
    
    !----------------------------------------
    write(*,*)   "Output_SST"
    
    open(fileid,file=trim(caseRoute)//"/"//gridNum//"/"//schemeNum//"/steady/"//trim(filename)//".dat")
    
    write(fileid,*) 'TITLE="euler solver on unstructured grid" '
    write(fileid,*) 'VARIABLES= "X" "Y"  "P" "Rou" "U"  "V" "muT" "K" "omg"' 
    write(fileid,*) "ZONE NODES=",nnodes," ELEMENTS=",ncells," ZONETYPE=FEQUADRILATERAL"
    write(fileid,*) "DATAPACKING=BLOCK"
    write(fileid,*) "VARLOCATION=([3-9]=CELLCENTERED)"
     
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
        write(fileid,*) muT(i)
    end do
    do i=1,ncells
        write(fileid,*)  U_SST(1,i)
    end do
    do i=1,ncells
        write(fileid,*)  U_SST(2,i)
    end do
    do i=1,ncells
        write(fileid,*) icell(:,i)
    end do
    close(fileid)
    !--------------------------------------------------------
    

end subroutine 
