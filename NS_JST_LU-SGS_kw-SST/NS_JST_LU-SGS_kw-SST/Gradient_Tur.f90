  subroutine Gradient_Tur
    implicit none
    integer::i,j
    integer::ncl,ncr
    
    Grad_UT = 0.0
    !the first iteration along the edge,calculate the gradient of u,v,T
    do i=1,nedges
      
        ncl=iedge(3,i)                      !atain the number of the left and the right cell
        ncr=iedge(4,i)
         
        if ( ncr .GT. 0 ) then 
            do j = 1,2
                Grad_UT(:,j,ncl) = Grad_UT(:,j,ncl) + UT_av(j,i)*vector(:,i)
                Grad_UT(:,j,ncr) = Grad_UT(:,j,ncr) - UT_av(j,i)*vector(:,i)
            end do
        else
             do j =1,2
                Grad_UT(:,j,ncl) = Grad_UT(:,j,ncl) + UT_av(j,i)*vector(:,i) 
             end do   
        end if   
    end do
    
    do j=1,2
        Grad_UT(1,j,:) = Grad_UT(1,j,:)/vol
        Grad_UT(2,j,:) = Grad_UT(2,j,:)/vol   
    end do
    
 end subroutine  
  
