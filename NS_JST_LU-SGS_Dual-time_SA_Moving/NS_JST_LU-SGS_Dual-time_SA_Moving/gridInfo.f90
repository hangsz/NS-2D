module gridInfo
    !--------------------------------------------------------------
        !purpose: read the grid infomation and calculate the values which are necessary for flow caluclation
    !--------------------------------------------------------------
    use controlPara    
    implicit none
    
    integer::nnodes,ncells,nedges 
    integer,allocatable:: iedge(:,:),icell(:,:)
    real(8),allocatable:: xy(:,:),vol(:)
    
    real(8),allocatable:: vector(:,:),ds(:)             
    real(8),allocatable:: cell_center(:,:)
    real(8),allocatable:: rL(:,:),rR(:,:)
    real(8),allocatable:: rij(:,:),lij(:),tij(:,:)
    
    real(8),allocatable::d(:)
    
    !---------------------------------------------------------------
    
               !variable specification  
               
    !nnodes     :   node numbers
	!ncells     :   cell numbers
	!nedges     :   edge numbers
	!iedge      :   matrix for edge
		            !iedge(1,i) : edge's start point 'a'
		            !iedge(2,i) : edge's end point 'b'
		            !iedge(3,i) : edge's left cell  (all positive)
		            !iedge(4,i) : positive for right cell
		                        !iedge(4,i) -1 for wall boundary
		                        !iedge(4,i) -2 for farfiled boundary
	!xy         :   cartesian coordinates  xy(1,i) for x  xy(2,i) for y
	!icell      :   triangle cell's three nodes index
	!vol        :   cell's volume(area in 2d)
    !vector     :   the normal vector of each edge
    !ds         :   the length of each edge
    !cell_center:   the (x,y) coordinates of each cell
    !rL,rR      :   the vectors points from the edge center to the left and right cell center respectively
    !rij,lij,tij:   the vector points from left cell center to the right; the magnitude of the vector; the unitized vector
    !d          :   the nearest distance of each cell to the wall  
    
    !--------------------------------------------------------------------
    
contains


subroutine grid

    call readGrid
    call gridData
    
    call distance
    
end subroutine

subroutine  readGrid
    !----------------------------------------------
         ! purpose : read the grid information
    !----------------------------------------------
         
    implicit none
	integer::i	
	
    integer::idGrid
    character(len = 200)::filename
    logical::alive
    
    !---------------------------------------------
    filename = trim(caseRoute)//"/"//gridNum//"/grid/naca0012.grd"
    
    inquire(file =filename, exist =alive)
    if ( .NOT. alive)  then
        write(*,*)  filename,'not exist.'
        stop
    end if
    
    open(newunit=idGrid,file=filename)
            
        read(idGrid,*)  nnodes,nedges,ncells
            
        ! allocate memory
        allocate(iedge(4,nedges))
        allocate(xy(2,nnodes))
        allocate(icell(4,ncells))
        allocate(vol(ncells))
            
        do i=1,nnodes
            read(idGrid,*)   xy(:,i)
        end do
            
        do i=1,nedges
            read(idGrid,*) iedge(:,i)
        end do
        
        do i=1,ncells
            read(idGrid,*) icell(:,i)
        end do
            
        do i=1,ncells
            read(idGrid,*) vol(i)
        end do
    close(idGrid)

end subroutine
 
subroutine  gridData
    !------------------------------------------------
        ! purpose: calculate related geometric variables
    !------------------------------------------------
         
    implicit none
    integer::i
     
    integer::start_node,end_node
    integer::ncl,ncr
     
    !------------------------------------------------
     
    allocate(vector(2,nedges))              
    allocate(ds(nedges))                    
    allocate(cell_center(2,ncells))
    allocate(rL(2,nedges))
    allocate(rR(2,nedges))
    allocate(rij(2,nedges))
    allocate(lij(nedges))
    allocate(tij(2,nedges))
     
     
    ! the normal vector and the length of every edge
    vector(1,:)= xy(2,iedge(2,:)) - xy(2,iedge(1,:))
    vector(2,:)=-xy(1,iedge(2,:)) + xy(1,iedge(1,:))
    ds(:)=sqrt( vector(1,:)**2+vector(2,:)**2 )
     
    cell_center = ( xy(:,icell(1,:) ) + xy(:,icell(2,:) ) + xy(:,icell(3,:) ) + xy(:,icell(4,:) ) ) /4.0
     
    do i=1,nedges
         
        start_node = iedge(1,i)
        end_node = iedge(2,i)
        ncl = iedge(3,i)
        ncr = iedge(4,i)
         
        select case(ncr)
             
        case(-2:-1)
            ! it's used to calculate the gradient of every edge
            rL(:,i) = 1.0/2*( xy(:,start_node) + xy(:,end_node) )  - cell_center(:,ncl)
            rij(:,i) = 2.0 * rL(:,i)
             
        case default
         
            rL(:,i) = 1.0/2*( xy(:,start_node) + xy(:,end_node) )  - cell_center(:,ncl)
            rR(:,i) = 1.0/2*( xy(:,start_node) + xy(:,end_node) )  - cell_center(:,ncr)
            rij(:,i) =  cell_center(:,ncr)-cell_center(:,ncl)   
             
        end select
                 
    end do
     
    lij= sqrt(rij(1,:)**2 + rij(2,:)**2 )
    tij(1,:) = rij(1,:) /lij
    tij(2,:) = rij(2,:) /lij
          
end subroutine
 
subroutine distance
    !---------------------------------------------------
         ! purpose: calculate the nearest distance of each cell to the wall
    !---------------------------------------------------
    
    implicit none
    integer::i,j
    integer::ncl,ncr
    integer::num = 0
    real(8)::d_temp
    real(8),allocatable::coor_b(:,:)
    
    !--------------------------------------------------
    
    allocate( d(ncells) )
    
    do i=1,nedges
        
        ncr = iedge(4,i)
        
        if(ncr .EQ. -1) num = num + 1
             
    end do
    
    allocate( coor_b(2,num) )
    j=0
    
    do i=1,nedges
        
        ncr = iedge(4,i)
        
        if(ncr .EQ. -1)  then
            
            j = j + 1
            coor_b(:,j) =  0.5*( xy(:,iedge(1,i) )  + xy(:,iedge(2,i) ) ) 
            
        end if
        
     end do
     
     do i =1,ncells
         
         do j=1,num
             
             if(j .EQ. 1)  then
                 
                 d(i)   = sqrt( ( cell_center(1,i) - coor_b(1,j) )**2 + ( cell_center(2,i) - coor_b(2,j) )**2 )
             else
                 
                 d_temp = sqrt( ( cell_center(1,i) - coor_b(1,j) )**2 + ( cell_center(2,i) - coor_b(2,j) )**2 )
                 if ( d_temp .LT. d(i) )    d(i)  = d_temp
             end if
           
         end do
     end do
     
end subroutine

 
end module