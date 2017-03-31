module moving
    use gridInfo
    implicit none
    
    real(8)::angularFrequency
    real(8),allocatable::xy0(:,:),vector0(:,:),rij0(:,:),rR0(:,:),rL0(:,:),tij0(:,:)
    
    !real(8),allocatable::U_Rot(:,:)
    !------------------------------------------------------
                ! variables specification
    ! xy0,vector0,...            :      stores the initial grid information before rotating.         
    ! U_Rot                      :      rotating speed of each node.
    !-----------------------------------------------------
    contains
    
subroutine allocateMemory_moving
    implicit none
    
    allocate( xy0(2,nnodes) )
    allocate( vector0(2,nedges) )
    allocate( rij0(2,nedges) ) 
    allocate( rR0(2,nedges) ) 
    allocate( rL0(2,nedges) )
    allocate( tij0(2,nedges) )
    
    !allocate( U_Rot(2,nnodes) )
    
end subroutine

subroutine  MovingInit
    implicit none
    ! store the initial grid information
    integer::i
    
    call allocateMemory_moving
    
    vector0 = vector
    rij0 =rij
    rR0 = rR
    rL0 = rL
    tij0 = tij
    
    do i=1,nnodes
        xy0(:,i) = xy(:,i)-rot_cen(:)    ! translate the grid node's coordinates to where the rotational center's coordinates are zero.
    end do 
    
    !-------------------------------------------------
    ! angular frequency
     angularFrequency =  2.0*kr*Ma_inf*sqrt(gamma*p_inf/rou_inf) !a_inf
    !-------------------------------------------------  
     
     
end subroutine

subroutine rotation( t,U_Rot,AOA ) 
    !------------------------------------------
       ! purpose:  renew all the geometric information with grid rotating. 
    !------------------------------------------
       
    implicit none 
    real(8),intent(in):: t
    real(8),intent(out)::U_Rot(:,:)
    real(8),intent(out)::AoA                     ! attact angle              
    real(8)::angular_rate            ! angular velocity 
    
    
    AoA = att + att_ampl* sin( angularFrequency * t )
    angular_rate = angularFrequency * ( att_ampl/180.0 *pi ) *cos( angularFrequency * t)
    
    ! all the geometric information has to be renewed with grid rotating
     
    xy(1,:) =   xy0(1,:) * cosd(AoA-att) + xy0(2,:) * sind(AoA-att) 
    xy(2,:) = - xy0(1,:) * sind(AoA-att) + xy0(2,:) * cosd(AoA-att) 
         
    U_Rot(1,:) =  angular_rate *  xy(2,:)   
    U_Rot(2,:) = -angular_rate *  xy(1,:)      
              
  
    vector(1,:) =  vector0(1,:) * cosd( AoA-att ) + vector0(2,:) * sind( AoA-att )
    vector(2,:) = -vector0(1,:) * sind( AoA-att ) + vector0(2,:) * cosd( AoA-att )   
    
    rij(1,:) =  rij0(1,:) * cosd( AoA-att ) + rij0(2,:) * sind( AoA-att )
    rij(2,:) = -rij0(1,:) * sind( AoA-att ) + rij0(2,:) * cosd( AoA-att )
    
    rR(1,:) =  rR0(1,:) * cosd( AoA-att ) + rR0(2,:) * sind( AoA-att )
    rR(2,:) = -rR0(1,:) * sind( AoA-att ) + rR0(2,:) * cosd( AoA-att )
    
    rL(1,:) =  rL0(1,:) * cosd( AoA-att ) + rL0(2,:) * sind( AoA-att )
    rL(2,:) = -rL0(1,:) * sind( AoA-att ) + rL0(2,:) * cosd( AoA-att )

    tij(1,:) =  tij0(1,:) * cosd( AoA-att ) + tij0(2,:) * sind( AoA-att )
    tij(2,:) = -tij0(1,:) * sind( AoA-att ) + tij0(2,:) * cosd( AoA-att )
              
end subroutine

       
end module
    