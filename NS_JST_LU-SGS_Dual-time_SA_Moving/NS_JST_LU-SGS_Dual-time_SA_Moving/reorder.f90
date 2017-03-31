subroutine reorder  
    !-----------------------------------------------
        ! purpose: reorder the cell number sequence
    !-----------------------------------------------
    
    implicit none
    integer:: i,j,k,kk
    integer::ncl,ncr
    integer::flag
    integer::cell_num,plane1
    
    integer::plane
    integer::member,member_last,member_next
    integer::cell_last,cell_next
    
    integer,allocatable::plane_cell(:)
    integer,allocatable::adjacent(:,:)
    integer,allocatable::assigned_last(:),assigned_next(:)
    
    character(len=200)::filename
    integer::idReo
    logical::alive1,alive2
    
    !----------------------------------------------  
    
    write(*,*)  "reorder"
    
    ! if stored file exist, read it directly
    filename = trim(caseRoute)//"/"//gridNum//"/"//schemeNum
    
    inquire(file = trim(filename)//"/LU-SGS/Assigned.dat",exist = alive1)
    inquire(file = trim(filename)//"/LU-SGS/Lower_upper.dat",exist = alive2)
    
    if( alive1 .AND. alive2 .AND. isNewGrid=="n" ) then
        open(newunit=idReo,file = trim(filename)//"/LU-SGS/Assigned.dat")
    
            do i=1,ncells
                 read(idReo,*)  assigned(i)
            end do
    
        close(idReo)
    
        open(newunit=idReo,file = trim(filename)//"/LU-SGS/Lower_Upper.dat") 
    
            do i =1,ncells
                read(idReo,*) lower(:,i)
            end do
            do i =1,ncells
                read(idReo,*) upper(:,i)
            end do
            
        close(idReo)  
        
        return  
     
    end if
    
    !-----------------------------------------------
    !no file exitst
    allocate( plane_cell(ncells) )                       ! the same plane 
    allocate( adjacent(4,ncells) )                        
    allocate( assigned_last(ncells) )
    allocate( assigned_next(ncells) )
    
    
    !-----------------------------------------------
    !find all adjacent cells of each cell
    adjacent = 0                                  ! it's set to be zero which can be useful to successive operation
    do i=1,nedges
        
        ncl = iedge(3,i)
        ncr = iedge(4,i)
        
        if( ncr.GT. 0) then                        ! only when the right cell exists is it necessary to continue
            do j=1,4                              ! three means that the cell is triangle, so it only has three edges
                if( adjacent(j,ncl) .EQ. 0 ) then  ! find the first number that equals zero,the cell number store here
                     adjacent(j,ncl) = ncr         ! left cell's adjacent cell is the right cell
                     exit                          ! if stored succesfully, then  exit the loop
                end if
            end do
            
            do j=1,4
                if( adjacent(j,ncr) .EQ. 0 ) then  
                     adjacent(j,ncr) = ncl
                     exit
                end if
            end do
           
        end if
    end do
        
    !----------------------------------------------
                         ! Reordering
    
    assigned = 0               !initialize
    assigned_last = 0
    assigned_next = 0
    
    !----------------------------------------------
    ! find the first cell on the corner or the wall and initialize all the arrays.
    do i=1,nedges
        ncl = iedge(3,i)
        ncr = iedge(4,i)
        if( ncr .LT. 0) exit
    end do
    
    plane_cell(ncl) = 1    ! the first cell member is 1 and plane also is 1
    
    assigned(1) = ncl      ! store the first cell number
    assigned_next(1) = ncl ! ncl is assigned this time ,it belongs to the assigned_next array.
    member_next = 1
    
    !----------------------------------
    plane = 1  
    member = 1
    
    do   
        if ( member .EQ. ncells)  exit   ! all the cell has been reordered, exit.
         
        member_last = member_next             !
        assigned_last = assigned_next         !
        plane = plane + 1
        
        member_next = 0
        assigned_next = 0
        
        do i=1,member_last                 ! the adjacent cell of the cell reordered previously are those to be reordered.
            
            cell_last = assigned_last(i)   ! the cell reordered last time
            
            do j=1,4 
                cell_next = adjacent(j,cell_last)       ! the cell is to be reordered
               
                if( cell_next .EQ. 0 ) then            ! this cell has no adjacent cell
                    exit                
                end if
                
                !---------------------
                ! judge wheather the cell being to be reordered has been reordered
                flag = 0
                do k=1,member 
                    if( cell_next .EQ. assigned(k)  )   then
                        flag = 1
                        exit
                    end if   
                end do
                
                if (flag .EQ. 1) cycle             ! reorder next adjacent cell
                
                !----------------------
                ! judge wheather the neighbers of the cell being to be reordered has been reordered to the same layer
                flag = 0
                
                do k=1,4                   ! all the neighbers of the cell being to reordered should be tested
                    if( adjacent(k,cell_next) .EQ. 0 ) exit  ! it states that all the adjacent cells have been tested
                   
                    do kk=1,member_next
                        if( adjacent(k,cell_next) .EQ. assigned_next(kk) ) then
                            flag = 1
                            exit
                        end if  
                    end do 
                    
                    if (flag .EQ. 1) exit
                    
                end do
             
                if (flag .EQ.1) cycle ! the neighber of the cell to be reordered has been reordered to the same layer.
                
                !----------------------------------------
                ! the cell has not been reordered
                member_next = member_next + 1
                member = member + 1   
                
                ! plane_cell(cell_next) = member_next
                plane_cell(cell_next) = plane
                
                assigned_next(member_next) = cell_next
                
                assigned(member)= cell_next
                
            end do   
         
        end do 
   
    end do
   
    !---------------------------------------------------
    ! upper and lower arrays store the edge  number respectively
    lower = 0
    upper = 0
    
    do i=1,nedges
        ncl = iedge(3,i)
        ncr = iedge(4,i)
        
        if (ncr .LT. 0)    cycle ! the edge has no right cell.
        
        plane1 = plane_cell(ncl)
        plane  = plane_cell(ncr)
        
        if( plane1 .LT. plane ) then                  ! the right cell in the left cell's upper array
            do k=1,4                                  ! the ncl is stored  where the value is zero
                if( upper(k,ncl) .EQ. 0 ) then
                    upper(k,ncl) = i                  ! store the edge number
                    exit
                end if
            end do 
            
            do k=1,4
                if( lower(k,ncr) .EQ. 0 ) then       ! same as the above
                    lower(k,ncr) = i
                    exit
                end if    
            end do
            
        else
            
            do k=1,4
                if( lower(k,ncl) .EQ. 0 ) then
                    lower(k,ncl) = i
                    exit
                end if      
            end do
            
            do k=1,4
                if( upper(k,ncr) .EQ. 0 ) then
                    upper(k,ncr) = i
                    exit
                end if    
            end do
            
        end if
        
    end do
  
    !-------------------------------------------------------
    !output to the file
    
    open(newunit=idReo,file = trim(filename)//"/LU-SGS/Assigned.dat")
    
        do i=1,ncells
            write(idReo,*)  assigned(i)
        end do
    
    close(idReo)
    
    open(newunit=idReo,file = trim(filename)//"/LU-SGS/Lower_Upper.dat")
    
        do i =1,ncells
            write(idReo,*) lower(:,i)
        end do
        do i =1,ncells
            write(idReo,*) upper(:,i)
        end do
    
    close(idReo)
    !------------------------------------------------------------
                             
end subroutine
