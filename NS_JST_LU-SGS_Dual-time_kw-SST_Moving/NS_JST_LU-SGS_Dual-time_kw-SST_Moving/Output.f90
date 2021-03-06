subroutine output(filename)
    !----------------------------------------------------
        ! purpose: out the flow properties to file
    !----------------------------------------------------
    implicit none 
    integer::i
    integer::idOut
    character(len=30)::filename
    
    !----------------------------------------------------
    !write(*,*)   "output"
    open(newunit = idOut,file=trim(caseRoute)//"/"//gridNum//"/"//schemeNum//"/flow/"//trim(filename)//".dat")
    
    write(idOut,*) 'TITLE="euler solver on unstructured grid" '
    write(idOut,*) 'VARIABLES= "X" "Y"  "P" "Rou" "U"  "V"  "nuT"' 
    write(idOut,*) "ZONE NODES=",nnodes," ELEMENTS=",ncells," ZONETYPE=FEQUADRILATERAL"
    write(idOut,*) "DATAPACKING=BLOCK"
    write(idOut,*) "VARLOCATION=([3-7]=CELLCENTERED)"
     
    do i=1,nnodes
        write(idOut,*) xy(1,i)
    end do
    do i=1,nnodes
        write(idOut,*) xy(2,i)
    end do
    
    do i=1,ncells
        write(idOut,*) U(5,i)
    end do
    do i=1,ncells
        write(idOut,*) U(1,i)
    end do
    do i=1,ncells
        write(idOut,*) U(2,i)
    end do
    do i=1,ncells
        write(idOut,*) U(3,i)
    end do
    
     do i=1,ncells
        write(idOut,*) muT(i)
    end do
    do i=1,ncells
        write(idOut,*) UT(1,i)
    end do
    do i=1,ncells
        write(idOut,*) UT(2,i)
    end do
    do i=1,ncells
        write(idOut,*) icell(:,i)
    end do
    close(idOut)
    
end subroutine 
