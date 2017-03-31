  subroutine gradient_tur 
    !----------------------------------------------------
        ! purpose: calculate the gradient of flow properties in the cell center. 
                 ! using Green's Theorem.
    !----------------------------------------------------
    
    implicit none
    integer::i,j
    integer::ncl,ncr
    
    real(8)::WT_av
  
    !---------------------------------------------------
    
    Grad_nuT = 0.0
    Grad_WT=0.0
  
    do i=1,nedges
      
        ncl=iedge(3,i)                   
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
  
