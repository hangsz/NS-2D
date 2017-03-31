  subroutine gradient                      
    !----------------------------------------------
        ! purpose: calcualate the gradient of the flow properties in the cell center.
                 ! using Green's Theorem.
    !----------------------------------------------
        
    implicit none
    integer::i,j
    integer::ncl,ncr
    
    Grad = 0.0
    
    do i=1,nedges
        
        ncl=iedge(3,i)                     
        ncr=iedge(4,i)
        
        if ( ncr .GT. 0 ) then 
            
            do j=1,6  
              Grad(:,j,ncl) = Grad(:,j,ncl) + U_av(j,i)*vector(:,i)
              Grad(:,j,ncr) = Grad(:,j,ncr) - U_av(j,i)*vector(:,i)
            end do 
            
        else
            
            do j=1,6  
                Grad(:,j,ncl) = Grad(:,j,ncl) + U_av(j,i)*vector(:,i) 
            end do 
            
        end if 
        
    end do
    
    do j=1,6
       Grad(1,j,:) = Grad(1,j,:)/vol
       Grad(2,j,:) = Grad(2,j,:)/vol
    end do
    
 end subroutine  
  
