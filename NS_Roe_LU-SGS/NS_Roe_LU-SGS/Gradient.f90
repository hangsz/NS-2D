  subroutine Gradient
    implicit none
    integer::i,j
    integer::ncl,ncr
    !real(8)::U_av(6)  !rou,u,v,w,p,T
    !allocate memory
    
    Grad = 0.0
    
    !the first iteration along the edge,calculate the gradient of u,v,T
    do i=1,nedges
        
        ncl=iedge(3,i)                      !atain the number of the left and the right cell
        ncr=iedge(4,i)
        if ( ncr .GT. 0 ) then 
            do j=1,6  
              !the gradient of u,v,T in every cell,left plus and right subtract
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
  
