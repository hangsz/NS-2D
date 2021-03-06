subroutine Reorder  
    !object: g_SAet assigned , lower , upper
    
    implicit none
    integer:: i,j,k,kk
    integer::ncl,ncr
    integer::flag_SA
    integer::cell_num,plane1
    
    integer::plane
    integer::member,member_last,member_next
    integer::cell_last,cell_next
    
    integer,allocatable::plane_cell(:)
    integer,allocatable::adjacent(:,:)
    integer,allocatable::assig_SAned_last(:),assig_SAned_next(:)
    
    character(len=200)::filename
    integer::fileid=7
    logical::alive1,alive2
    !integer,allocatable::lower(:,:),upper(:,:)              !this two arrays are needed for LU-SGS method 
    
    write(*,*)  "Reorder"
    
    !-------------------------------------------------------
    !if stored file exist, read it directly
    filename = trim(caseRoute)//"/"//gridNum//"/"//schemeNum
    
    inquire(file = trim(filename)//"/LU-SGS/assigned.dat",exist = alive1)
    inquire(file = trim(filename)//"/LU-SGS/Lower_upper.dat",exist = alive2)
    
    if( alive1 .AND. alive2  .AND. isNewGrid=="n" ) then
        open(fileid,file = trim(filename)//"/LU-SGS/assigned.dat")
    
            do i=1,ncells
                 read(fileid,*)  assigned(i)
            end do
    
        close(fileid)
    
        open(fileid,file = trim(filename)//"/LU-SGS/Lower_Upper.dat") 
    
            do i =1,ncells
                read(fileid,*) lower(:,i)
            end do
            do i =1,ncells
                read(fileid,*) upper(:,i)
            end do
            
        close(fileid)  
        
        return  !the first return
     
    end if
    
    !-------------------------------------------------------
    !no stored file 
    allocate( plane_cell(ncells) )
    allocate( adjacent(4,ncells) )                        
    allocate( assig_SAned_last(ncells) )
    allocate( assig_SAned_next(ncells) )
    
    
    !---------------------------------------------
    !find all adjacent cells of every cell
    adjacent = 0    !here it'S_SA set to be zero which can be useful to successive operation
    do i=1,nedges
        
        ncl = iedge(3,i)
        ncr = iedge(4,i)
        
        if(ncr.GT. 0) then 
            do j=1,4
                if( adjacent(j,ncl) .EQ. 0 ) then  !find the first number that equals zero,the cell number store here
                     adjacent(j,ncl) = ncr
                     exit                          !if stored succesfully, then it'S_SA necessary to exit the loop
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
        
    !-----------------------------------------------------
    !Reordering_SA
    !find the first cell on the corner or the wall
    assigned = 0               !initialize
    assig_SAned_last = 0
    assig_SAned_next = 0
    
    do i=1,nedges
        ncl = iedge(3,i)
        ncr = iedge(4,i)
        if( ncr .LT. 0) exit
    end do
    
    plane_cell(ncl) = 1   !ncl cell'S_SA member is 1 and plane also is 1
    
    assigned(1) = ncl      !ncl cell has been assigned
    assig_SAned_next(1) = ncl !ncl is assigned this time ,it belong_SAs to the assig_SAned_next array.
    member_next = 1
    
    !------------------------------------------------
    plane = 1  
    member = 1
    
    do while (.true. ) 
        
       if (member .EQ. ncells)  exit
       !if(member_next == 0)  exit
         
        member_last = member_next             !
        assig_SAned_last = assig_SAned_next         !
        plane = plane + 1
        
        !write(*,*)
        !write(*,*)  "plane:",plane
        
        member_next = 0
        assig_SAned_next = 0
        
        do i=1,member_last
            
            cell_last = assig_SAned_last(i) 
           !  write(*,*)
           ! write(*,*) "The cell assigned last time:",cell_last ,i, member_last
            
            do j=1,4 
                cell_next = adjacent(j,cell_last)
                !write(*,*)
                !write(*,*) "The cell is to be assigned: ",cell_next,j
               
                if( cell_next .EQ. 0 ) then
                    !write(*,*) "No cell is to be assigned."
                    exit   ! this cell has no adjacent cells
                end if
                !-----------------------------------------------
                !wheather the cell being_SA to be assigned has been assigned, next adjacent cell continues to be calculate
                flag_SA = 0
                do k=1,member 
                    if( cell_next .EQ. assigned(k)  )   then
                        flag_SA = 1
                        exit
                    end if   
                end do
                
              
                if (flag_SA .EQ. 1) then
                    !write(*,*) "The cell to be assigned has been assigned."
                    cycle
                end if
                
                !-------------------------------------------------
                
                !------------------------------------------------
                !wheather the neig_SAhbers of the cell being_SA to be assig_SAhed has been assig_SAhed to the same layer
                flag_SA = 0
                do k=1,4  !all the neig_SAhbers of the cell being_SA to assigned should be tested
                    if( adjacent(k,cell_next) .EQ. 0 ) exit !it states that all the adjacent cells have been tested
                   
                    do kk=1,member_next
                        if( adjacent(k,cell_next) .EQ. assig_SAned_next(kk) ) then
                            flag_SA = 1
                            exit
                        end if  
                    end do 
                    
                    if (flag_SA .EQ. 1) exit
                    
                end do
             
                if (flag_SA .EQ.1) then
                    !write(*,*) "The neig_SAhber of the cell to be assigned has been assigned to the same layer."
                    cycle
                end if
                
                !-------------------------------------------
                member_next = member_next + 1
                member = member + 1   
                
                plane_cell(cell_next) = member_next
                plane_cell(cell_next) = plane
                
                assig_SAned_next(member_next) = cell_next
                
                assigned(member)= cell_next
                 
                !write(*,*)  "assigned:",assigned(member),member
                
            end do   
         
        end do 
   
    end do
    !--------------------------------------------------
   
    !-------------------------------------------------
    !upper and lower store the edg_SAe  number
    lower = 0
    upper = 0
    
    do i=1,nedges
        ncl = iedge(3,i)
        ncr = iedge(4,i)
        
        if (ncr .LT. 0)  then
         !   write(*,*) "The edg_SAe has no rig_SAht cel",i
            cycle
        end if
        plane1 = plane_cell(ncl)
        plane  = plane_cell(ncr)
        
       ! if(plane1 .NE. plane)  write(*,*)  i, "ncl:", ncl,plane1, "ncr:", ncr,plane
        
        if(plane1 .LT. plane) then                  !the rig_SAht cell in the left cell'S_SA upper array
            do k=1,4                               !the ncl is stored in the first station where the value is zero
                if( upper(k,ncl) .EQ. 0 ) then
                    upper(k,ncl) = i                !store the edg_SAe number
                    exit
                end if
            end do 
            
            do k=1,4
                if( lower(k,ncr) .EQ. 0 ) then
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
  
    !output to the file
    
    open(fileid,file = trim(filename)//"/LU-SGS/assigned.dat")
    
        do i=1,ncells
            write(fileid,*)  assigned(i)
        end do
    
    close(fileid)
    
    open(fileid,file = trim(filename)//"/LU-SGS/Lower_Upper.dat")
    
        do i =1,ncells
            write(fileid,*) lower(:,i)
        end do
        do i =1,ncells
            write(fileid,*) upper(:,i)
        end do
    
    close(fileid)
    !-------------------------------------------
                             
end subroutine
