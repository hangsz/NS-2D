  subroutine Gradient_Tur
    implicit none
    integer::i,j
    integer::ncl,ncr
    
    real(8)::WT_av
    !allocate memory
    
    Grad_nuT = 0.0
    Grad_WT=0.0
    !the first iteration along the edge,calculate the gradient of u,v,T
    do i=1,nedges
      
        ncl=iedge(3,i)                      !atain the number of the left and the right cell
        ncr=iedge(4,i)
        
        WT_av = U_av(1,i) *nuT_av(i)
        if ( ncr .GT. 0 ) then 
            Grad_nuT(:,ncl) = Grad_nuT(:,ncl) + nuT_av(i)*vector(:,i)
            Grad_nuT(:,ncr) = Grad_nuT(:,ncr) - nuT_av(i)*vector(:,i)
            Grad_WT(:,ncl) = Grad_WT(:,ncl) + WT_av*vector(:,i)
            Grad_WT(:,ncr) = Grad_WT(:,ncr) - WT_av*vector(:,i)
        else
            Grad_nuT(:,ncl) = Grad_nuT(:,ncl) + nuT_av(i)*vector(:,i) 
            Grad_WT(:,ncl) = Grad_WT(:,ncl) + WT_av*vector(:,i) 
        end if 
    end do
    
       Grad_nuT(1,:) = Grad_nuT(1,:)/vol
       Grad_nuT(2,:) = Grad_nuT(2,:)/vol
       Grad_WT(1,:) = Grad_WT(1,:)/vol
       Grad_WT(2,:) = Grad_WT(2,:)/vol
       
 end subroutine  
  
