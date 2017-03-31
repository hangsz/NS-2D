module Grid_info
    
    use Control_para
    
    implicit none
    integer::nnodes,ncells,nedges 
    integer,allocatable:: iedge(:,:),icell(:,:)
    real(8),allocatable:: xy(:,:),vol(:)
    real(8),allocatable:: vector(:,:),ds(:)             
    real(8),allocatable:: cell_center(:,:)
    real(8),allocatable:: rL(:,:),rR(:,:)
    real(8),allocatable:: rij(:,:),lij(:),tij(:,:)
    
    real(8),allocatable::d(:)
    
contains


subroutine Grid
    call Read_grid
    call Grid_data
    
    call Distance
    
end subroutine

subroutine  Read_grid
    implicit none
	integer::i	
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

    open(10,file= trim(caseRoute)//"/"//gridNum//"/grid/naca0012.grd")
      
        read(10,*)  nnodes,nedges,ncells
            
        !allocate memory
        allocate(iedge(4,nedges))
        allocate(xy(2,nnodes))
        allocate(icell(4,ncells))
        allocate(vol(ncells))
            
        do i=1,nnodes
            read(10,*)   xy(:,i)
        end do
            
        do i=1,nedges
            read(10,*) iedge(:,i)
        end do
        
        do i=1,ncells
            read(10,*) icell(:,i)
        end do
            
        do i=1,ncells
            read(10,*) vol(i)
        end do
    close(10)

 end subroutine
 
 subroutine  Grid_data
     implicit none
     integer::i
     
     integer::start_node,end_node,ncl,ncr
     !------------------------------------------------------
     !allocate memory
     allocate(vector(2,nedges))               !the normal vector  of every edge
     allocate(ds(nedges))                     !the  length of every edge
     allocate(cell_center(2,ncells))
     allocate(rL(2,nedges))
     allocate(rR(2,nedges))
     allocate(rij(2,nedges))
     allocate(lij(nedges))
     allocate(tij(2,nedges))
     
     
     !the normal vector and the length of every edge, clockwise
     vector(1,:)= xy(2,iedge(2,:)) - xy(2,iedge(1,:))         ! dy
     vector(2,:)=-xy(1,iedge(2,:)) + xy(1,iedge(1,:))         ! -dx
     ds(:)=sqrt( vector(1,:)**2+vector(2,:)**2 )
     !-------------------------------------------------------
     
     cell_center = ( xy(:,icell(1,:) ) + xy(:,icell(2,:) ) + xy(:,icell(3,:) ) + xy(:,icell(4,:) ) ) /4.0
     
     do i=1,nedges
         start_node = iedge(1,i)
         end_node = iedge(2,i)
         ncl = iedge(3,i)
         ncr = iedge(4,i)
         select case(ncr)
         case(-2:-1)
             ! rL;the left cell center --> edge center
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
 
 subroutine Distance 
    !-----------------------------------------------
        ! purpose: calculate the distance of each cell to the wall
    !-----------------------------------------------
    implicit none
    integer::i,j
    integer::ncl,ncr
    integer::num = 0
    
    
    real(8)::d_temp
    real(8),allocatable::coor_b(:,:)
    
    allocate( d(ncells) )
    do i=1,nedges
        ncr = iedge(4,i)
        
        if(ncr .EQ. -1)  then
            num = num + 1
        end if
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
         d(i)   = sqrt( ( cell_center(1,i) - coor_b(1,1) )**2 + ( cell_center(2,i) - coor_b(2,1) )**2 ) 
         do j=2,num

            d_temp = sqrt( ( cell_center(1,i) - coor_b(1,j) )**2 + ( cell_center(2,i) - coor_b(2,j) )**2 )
            if ( d_temp .LT. d(i) )  then
                d(i)  = d_temp
            end if
                 
         end do
     end do
     
end subroutine

 
 end module