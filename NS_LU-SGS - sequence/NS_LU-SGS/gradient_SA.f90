  subroutine gradient_SA( U_av)
    implicit none
    real(8),intent(in)::U_av(:,:)
    integer::i,j
    integer::ncl,ncr
    
    real(8)::WT_av
    !allocate memory
    
    Grad_U_SA = 0.0
    Grad_W_SA=0.0
    !the first iteration along the edge,calculate the gradient of u,v,T
    do i=1,nedges
      
        ncl=iedge(3,i)                      !atain the number of the left and the right cell
        ncr=iedge(4,i)
        
        WT_av = U_av(1,i) *U_av_SA(i)
        if ( ncr .GT. 0 ) then 
            Grad_U_SA(:,ncl) = Grad_U_SA(:,ncl) + U_av_SA(i)*vector(:,i)
            Grad_U_SA(:,ncr) = Grad_U_SA(:,ncr) - U_av_SA(i)*vector(:,i)
            Grad_W_SA(:,ncl) = Grad_W_SA(:,ncl) + WT_av*vector(:,i)
            Grad_W_SA(:,ncr) = Grad_W_SA(:,ncr) - WT_av*vector(:,i)
        else
            Grad_U_SA(:,ncl) = Grad_U_SA(:,ncl) + U_av_SA(i)*vector(:,i) 
            Grad_W_SA(:,ncl) = Grad_W_SA(:,ncl) + WT_av*vector(:,i) 
        end if 
    end do
    
       Grad_U_SA(1,:) = Grad_U_SA(1,:)/vol
       Grad_U_SA(2,:) = Grad_U_SA(2,:)/vol
       Grad_W_SA(1,:) = Grad_W_SA(1,:)/vol
       Grad_W_SA(2,:) = Grad_W_SA(2,:)/vol
       
 end subroutine  
  
