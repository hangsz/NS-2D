subroutine output(filename)
    implicit none 
    integer::i
    integer,parameter::fileid=7
    character(len=50)::filename
    !----------------------------------------------------------
    write(*,*)   "Output"
    
    if( moveNum==0 )  then
        
        open(fileid,file=trim(caseRoute)//"/"//trim(gridNum)//"/"//trim(fileNum)//"/steady/"//trim(filename)//".dat")
    else
        open(fileid,file=trim(caseRoute)//"/"//trim(gridNum)//"/"//trim(fileNum)//"/flow/"//trim(filename)//".dat")
    end if
    
    select case( turModelNum) 
        
    case( 0 )  
        
        write(fileid,*) 'TITLE="Flow properties" '
        write(fileid,*) 'VARIABLES= "X" "Y"  "P" "Rou" "U"  "V" ' 
        write(fileid,*) "ZONE NODES=",nnodes," ELEMENTS=",ncells," ZONETYPE=FEQUADRILATERAL"
        write(fileid,*) "DATAPACKING=BLOCK"
        write(fileid,*) "VARLOCATION=([3-6]=CELLCENTERED)"
     
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
            write(fileid,*) icell(:,i)
        end do
        
    case( 1 )
        
        write(fileid,*) 'TITLE="Flow properties" '
        write(fileid,*) 'VARIABLES= "X" "Y"  "P" "Rou" "U"  "V" "muT" "nu"' 
        write(fileid,*) "ZONE NODES=",nnodes," ELEMENTS=",ncells," ZONETYPE=FEQUADRILATERAL"
        write(fileid,*) "DATAPACKING=BLOCK"
        write(fileid,*) "VARLOCATION=([3-8]=CELLCENTERED)"
     
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
            write(fileid,*)  muT(i)
        end do
        
        do i=1,ncells
            write(fileid,*)  U_SA(i)
        end do
        do i=1,ncells
            write(fileid,*) icell(:,i)
        end do
        
    case( 2 )
        
        write(fileid,*) 'TITLE="Flow properties" '
        write(fileid,*) 'VARIABLES= "X" "Y"  "P" "Rou" "U"  "V" "mu" "k" "omg" ' 
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
            write(fileid,*)  muT(i)
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
        
    case default
    
    end select
    
    
    close(fileid)
  
end subroutine 
