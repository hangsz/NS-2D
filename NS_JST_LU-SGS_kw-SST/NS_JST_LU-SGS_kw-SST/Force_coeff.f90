subroutine Force_coeff(filename)
    implicit none
    integer::num_u,num_l,num_b
    integer,allocatable::,wall(:),wall_u(:),wall_l(:)
    real(8)::Cp,Cn,Cm
    real(8)::temp
    integer,parameter::fileid=7
  
    integer::i,j,k
    logical::alive1,alive2
    character(len =30)::filename
    
    inquire(file="Wall_u.dat",exist = alive1 )
    inquire(file="Wall_l.dat",exist = alive2 )
    
    if( alive1 .AND. alive2 )  then
        
        open(fileid,file = "Wall_u.dat")
            
        read(fileid,*)  num_u
        
        allocate( wall_u(num_u) )  
        
        do i=1,num_u
                read(fileid,*)  wall_u(i)
         end do
    
        close(fileid)
        
        open(fileid,file = "wall_l.dat")
            
        read(fileid,*)  num_l
        
        allocate( wall_l(num_l) )  
        do i=1,num_l
                read(fileid,*)  wall_l(i)
         end do
    
        close(fileid)
        
    else
        !***************************
        !get the total numbers of cells on the wall
        num_b=0
        
        do i =1,nedges
            if(iedge(4,i) .EQ. -1) then 
                num_b = num_b + 1    
            end if  
        end do
     
        !get each cell number 
        allocate( wall(num_b) )  
        j=1
        do i =1,nedges 
            if( iedge(4,i) .EQ. -1 ) then 
                wall(j) = i
                j=j+1
            end if
        end do
        !***************************
        !divide the cells intoto the upper and lower array ,according to the y cooridnates
        num_u=0
        num_l=0
        do i=1,num_b
            ncl=iedge(3,wall(i) )
            if( cell_center(2, ncl ) .GT. 0.0 )  then
          
                num_u=num_u+1
            else
                num_l=num_l+1
            end if
        end do
        
        allocate( wall_u(num_u ) )
        allocate( wall_l(num_l ) )
        !divide the cells intoto the upper and lower array ,according to the y cooridnates
        j=1
        k=1
        do i=1,num_b
            ncl=iedge(3,wall(i) )
            if( cell_center(2, ncl ) .GT. 0.0 )  then
                wall_u(j)=wall(i)
                j=j+1
            else
                wall_l(k)=wall(i)
                k=k+1
            end if
        end do
     
       !**********************
        !output the numberof cells  on the wall
        open(fileid,file = "wall_u.dat")
            
        write(fileid,*)  num_u
        
         do i=1,num_u
                write(fileid,*)  wall_u(i)
         end do
    
        close(fileid)
        
        open(fileid,file = "wall_l.dat")
            
        write(fileid,*)  num_l
        
         do i=1,num_l
                write(fileid,*) wall_l(i)
         end do
    
        close(fileid)
        
    end if
   !*************************************************************************
    Cl=0.0
    Cm=0.0
    
    filename1 = trim(filename)//"Cn_Cm_att.dat"

    do i=1,num_u 
        ncl=iedge(3,wall_u(i) )
        Cp = ( U(5, ncl)-p_inf ) / ( 1.0/2*rou_inf*( u_inf**2 + v_inf**2 )  )
        
        
        Cn = Cl - Cp*ds( wall_u(i) ) *  cos( atan( vector(1,wall_u(i) ) / vector(2,wall_u(i) )  )  )   
       
        Cm = Cm + ( - Cp*ds( wall_u(i) ) *  cos( atan( vector(2,wall_u(i) ) / vector(1,wall_u(i) )  )  ) * ( cell_center(1,ncl) - rot_cen(1) )   &
               &+ ( - Cp*ds( wall_u(i) ) *  sin( atan( vector(2,wall_u(i) ) / vector(1,wall_u(i) )  )  ) * ( cell_center(2,ncl) - rot_cen(2) )
        
    end do
    
    do i=1,num_l
        ncl=iedge(3,wall_l(i) )
        Cp = ( U(5, ncl)-p_inf ) / ( 1.0/2*rou_inf*(u_inf**2 + v_inf**2 )  )
        
        Cn = Cl - Cp*ds( wall_u(i) ) *  cos( atan( vector(1,wall_u(i) ) / vector(2,wall_u(i) )  )  )   
       
        Cm = Cm +  Cp*ds( wall_u(i) ) *  cos( atan( vector(2,wall_u(i) ) / vector(1,wall_u(i) )  )   * ( cell_center(1,ncl) - rot_cen(1) )   &
               &+  Cp*ds( wall_u(i) ) *  sin( atan( vector(2,wall_u(i) ) / vector(1,wall_u(i) )  )   * ( cell_center(2,ncl) - rot_cen(2) )
       
    end do
  
    close(fileid)
    
    filename1 = trim(filename) //"-Cp_l.dat"
    open(fileid,file=filename1)
   
    
    
end subroutine

            
    
        
    
    
    
    
    
    
    