  subroUTine gradient_SST
    implicit none
    integer::i,j
    integer::ncl,ncr
    
    Grad_SST = 0.0
    !the first iteration along the edge,calculate the gradient of u,v,T
    do i=1,nedges
      
        ncl=iedge(3,i)                      !atain the number of the left and the right cell
        ncr=iedge(4,i)
         
        if ( ncr .GT. 0 ) then 
            do j = 1,2
                Grad_SST(:,j,ncl) = Grad_SST(:,j,ncl) + U_av_SST(j,i)*vector(:,i)
                Grad_SST(:,j,ncr) = Grad_SST(:,j,ncr) - U_av_SST(j,i)*vector(:,i)
            end do
        else
             do j =1,2
                Grad_SST(:,j,ncl) = Grad_SST(:,j,ncl) + U_av_SST(j,i)*vector(:,i) 
             end do   
        end if   
    end do
    
    do j=1,2
        Grad_SST(1,j,:) = Grad_SST(1,j,:)/vol
        Grad_SST(2,j,:) = Grad_SST(2,j,:)/vol   
    end do
    
 end subroUTine  
  
