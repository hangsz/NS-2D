subroutine Rotation(AoA,angular_rate)  !  the next step's node speed and coordinates
    implicit none 
    real(8)::AoA   !angle of attack
    real(8)::angular_rate
     
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
    !
    tij(1,:) =  tij0(1,:) * cosd( AoA-att ) + tij0(2,:) * sind( AoA-att )
    tij(2,:) = -tij0(1,:) * sind( AoA-att ) + tij0(2,:) * cosd( AoA-att )
              
end subroutine

    