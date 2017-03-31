subroutine Cp_cal(filename)
    implicit none
    integer::num_u,num_l,num_b
    integer,allocatable::cell_u(:),cell_l(:),cell_b(:)
    real(8)::Cp
    integer::temp
    integer,parameter::fileid=7
  
    integer::i,j,k
    logical::alive1,alive2
    character(len =30)::filename,filename1
    
    
    
    inquire(file="Cell_u.dat",exist = alive1 )
    inquire(file="Cell_l.dat",exist = alive2 )
    
    if( alive1 .AND. alive2 )  then
        
        open(fileid,file = "Cell_u.dat")
            
        read(fileid,*)  num_u
        
        allocate( cell_u(num_u) )  
        
        do i=1,num_u
                read(fileid,*)  cell_u(i)
         end do
    
        close(fileid)
        
        open(fileid,file = "Cell_l.dat")
            
        read(fileid,*)  num_l
        
        allocate( cell_l(num_l) )  
        
        do i=1,num_l
                read(fileid,*)  cell_l(i)
         end do
    
        close(fileid)
        
    else
        !----------------------------------------------
        !get the total numbers of cells on the wall
        num_b=0
        
        do i =1,nedges
            if(iedge(4,i) .EQ. -1) then 
                num_b = num_b + 1    
            end if  
        end do
     
        !get each cell number 
        allocate( cell_b(num_b) )  
        
        j=1
        do i =1,nedges 
            if( iedge(4,i) .EQ. -1 ) then 
                cell_b(j) = iedge(3,i) 
                !write(*,*) j,cell_b(j) 
                j=j+1
            end if
        end do
        !-----------------------------------
        !divide the cells intoto the upper and lower array ,according to the y cooridnates
        num_u=0
        num_l=0
        do i=1,num_b
            if( cell_center(2, cell_b(i) ) .GT. 0.0 )  then
          
                num_u=num_u+1
            else
                num_l=num_l+1
            end if
        end do
        
        
        allocate( cell_u(num_u ) )
        allocate( cell_l(num_l ) )
        !divide the cells intoto the upper and lower array ,according to the y cooridnates
        j=1
        k=1
        do i=1,num_b
            if( cell_center(2, cell_b(i) ) .GT. 0.0 )  then
                cell_u(j)=cell_b(i)
                j=j+1
            else
                cell_l(k)=cell_b(i)
                k=k+1
            end if
        end do
        !-----------------------------------------
        !sort the upper and lower cell according to the x coorinates
        do j=1,num_u-1
            do i=num_u-1,j,-1
            
                if( cell_center(1,cell_u(i) ) .GT. cell_center(1,cell_u(i+1) ) )  then
                      temp=cell_u(i)
                      cell_u(i) = cell_u(i+1)
                      cell_u(i+1) = temp
                end if
            end do
        end do
        do j=1,num_l-1
            do i=num_l-1,j,-1
            
                if( cell_center(1,cell_l(i) ) .GT. cell_center(1,cell_l(i+1) ) )  then
                      temp=cell_l(i)
                      cell_l(i) = cell_l(i+1)
                      cell_l(i+1) = temp
                end if
            end do
        end do
        
              
       !------------------------------------
        !output the numberof cells  on the wall
        open(fileid,file = "Cell_u.dat")
            
        write(fileid,*)  num_u
        
         do i=1,num_u
                write(fileid,*)  cell_u(i)
         end do
    
        close(fileid)
        
        open(fileid,file = "Cell_l.dat")
            
        write(fileid,*)  num_l
        
         do i=1,num_l
                write(fileid,*)  cell_l(i)
         end do
    
        close(fileid)
        
    end if
   !----------------------------------------------------
    
    filename1 = trim(filename)//"-Cp_u.dat"
    open(fileid,file=filename1)
    do i=1,num_u 
        Cp = ( U(5, cell_u(i) )-p_inf ) / ( 1.0/2*U(1, cell_u(i) )*( U(2, cell_u(i) )**2 + U(3, cell_u(i))**2 )  )
        write(fileid,*)  cell_center(:,cell_u(i) ),Cp
    end do
    
    close(fileid)
    
    filename1 = trim(filename) //"-Cp_l.dat"
    open(fileid,file=filename1)
    do i=1,num_l
        Cp = ( U(5, cell_l(i) )-p_inf ) / ( 1.0/2*U(1, cell_l(i) )*( U(2, cell_l(i) )**2 + U(3, cell_l(i) )**2 )  )
        write(fileid,*)  cell_center(:,cell_l(i) ),Cp
    end do
    close(fileid)
    
    
end subroutine

            
    
        
    
    
    
    
    
    
    