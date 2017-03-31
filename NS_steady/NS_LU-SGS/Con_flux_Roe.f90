subroutine Con_flux_Roe   !convective flux adopte Roe scheme
    implicit none
    integer::i,j
    integer::ncl,ncr
    
    real(8)::epsi2L,epsi2R,dhL,dhR
    real(8)::d1maxL(5),d1minL(5),d1maxR(5),d1minR(5)
    real(8)::UL(5),UR(5)  !rou,u,v,w,p
    
    real(8)::Vn                 ! the contravariant velocity- the velocity normal to the surface element
    real(8)::d2L(5),d2R(5)
    real(8)::psiL(5),psiR(5)
    real(8)::nx,ny
   
    real(8)::HL,HR,VnL,VnR
    real(8)::rou_Roe,u_Roe,v_Roe,H_Roe,q2_Roe,c_Roe,Vn_Roe
    real(8)::FcL(5),FcR(5),dF1(5),dF234(5),dF5(5)
    
    real(8)::delta  !HarT_SAen'S_SA entropy correction
    real(8)::lamda_c1,lamda_c234,lamda_c5  !|V_Roe - c_Roe|,|V_Roe|,|V_Roe + c_Roe|
    
    real(8),allocatable::Umax(:,:),Umin(:,:),psi(:,:)   

    !Venkatakrishnan'S_SA limiter  for Roe scheme
    allocate(Umax(5,ncells))  
    allocate(Umin(5,ncells))
    allocate(psi(5,ncells))
    
    !---------------------------------------------------------------------------------
    !calculathe the Umax Umin
    Umax= U
    Umin= U
    
    do i=1,nedges
        
        ncl=iedge(3,i)
        ncr=iedge(4,i)
        
        if ( ncr .GT. 0)  then 
          !the maximun/minimum value is among_SA the last Umax and the current left and rig_SAht valus
        
            do j=1,5
              
                Umax(j,ncl) = max( Umax(j,ncl),U(j,ncl),U(j,ncr) )
                Umin(j,ncl) = min( Umin(j,ncl),U(j,ncl),U(j,ncr) )
          
                Umax(j,ncr) = max( Umax(j,ncr),U(j,ncl),U(j,ncr) )
                Umin(j,ncr) = min( Umin(j,ncr),U(j,ncl),U(j,ncr) ) 
            end do
           
        else 
            do j=1,5  
                Umax(j,ncl) = max( Umax(j,ncl),U(j,ncl) )
                Umin(j,ncl) = min( Umin(j,ncl),U(j,ncl) )
            end do
              
        end if  
        
    
    end do
    
    !-------------------------------------------------------------------------------
   
    !Venkatakrishnan'S_SA limiter
    !notice : all quantities must be non-dimensional
    
    psi = 1.0
     
    do i=1,nedges
        
        ncl=iedge(3,i)
        ncr=iedge(4,i)
     
        if( ncr .GT. 0 ) then
            d1maxL =  Umax(:,ncl) - U(:,ncl) 
            d1minL =  Umin(:,ncl) - U(:,ncl)
            
            d1maxR =  Umax(:,ncr) - U(:,ncr)
            d1minR =  Umin(:,ncr) - U(:,ncr)
         
            dhL =sqrt( vol(ncl) )
            dhR =sqrt( vol(ncr) )
            
            epsi2L = ( K*dhL )**3
            epsi2R = ( K*dhR )**3
            
            do j=1,5   
               !psiL
                d2L(j) = dot_product( Grad(:,j,ncl),rL(:,i) ) 
                d2R(j) = dot_product( Grad(:,j,ncr),rR(:,i) )
                
                d2L(j) =sign( abs(d2L(j)) + 10.0E-12 , d2L(j) )
                d2R(j) =sign( abs(d2R(j)) + 10.0E-12 , d2R(j) )
                
                if ( d2L(j) .GT. 0.0 ) then
                    psiL(j) = 1.0/d2L(j)* ( (d1maxL(j)**2 + epsi2L)*d2L(j) + 2.0 * d2L(j)**2 * d1maxL(j) ) / ( d1maxL(j)**2 + 2.0*d2L(j)**2 + d1maxL(j)*d2L(j) + epsi2L )
                else if( d2L(j) .LT. 0.0 ) then
                    psiL(j) = 1.0/d2L(j)* ( (d1minL(j)**2 + epsi2L)*d2L(j) + 2.0 * d2L(j)**2 * d1minL(j) ) / ( d1minL(j)**2 + 2.0*d2L(j)**2 + d1minL(j)*d2L(j) + epsi2L )   
                else
                    psiL(j) = 1.0
                end if
                
                !psiR
                if ( d2R(j) .GT. 0.0 ) then
                    psiR(j) = 1.0/d2R(j)* ( (d1maxR(j)**2 + epsi2R)*d2R(j) + 2.0 * d2R(j)**2 * d1maxR(j) ) / ( d1maxR(j)**2 + 2.0*d2R(j)**2 + d1maxR(j)*d2R(j) + epsi2R )
                else if( d2R(j) .LT. 0.0 ) then
                    psiR(j) = 1.0/d2R(j)* ( (d1minR(j)**2 + epsi2R)*d2R(j) + 2.0 * d2R(j)**2 * d1minR(j) ) / ( d1minR(j)**2 + 2.0*d2R(j)**2 + d1minR(j)*d2R(j) + epsi2R )   
                else 
                    psiR(j) = 1.0
                end if
            
                psi(j,ncl)  = min( psi(j,ncl),psiL(j),psiR(j) )
                psi(j,ncr)  = min( psi(j,ncr),psiL(j),psiR(j) )
                
               ! if( psi(j,ncl) < 0.0  .OR. psi(j,ncl) >1.0 ) then
               !     stop "1"
               ! else
               !     write(*,*)  psi(j,ncl)
               ! end if
               ! 
               ! if( psi(j,ncr) < 0.0  .OR. psi(j,ncr) >1.0 ) then
               !     stop "2"
               ! else
               !     write(*,*)  psi(j,ncr)
               !end if
               
            end do 
    
        end if
        
    end do
    
    !-------------------------------------------------------------------------------
    !calculathe  Fc
    
    Fc=0.0
    
    do i=1,nedges
        ncl=iedge(3,i)
        ncr=iedge(4,i)
       
        select case(ncr)
        case( -2:-1 ) 
            
            Vn = dot_product( U_av(2:3,i),vector(:,i) )
            
            Fc(1,ncl) = Fc(1,ncl) + U_av(1,i)*Vn
            Fc(2,ncl) = Fc(2,ncl) + U_av(1,i)*U_av(2,i)*Vn + U_av(5,i)*vector(1,i)
            Fc(3,ncl) = Fc(3,ncl) + U_av(1,i)*U_av(3,i)*Vn + U_av(5,i)*vector(2,i)
            Fc(5,ncl) = Fc(5,ncl) + U_av(1,i)*Vn*( U_av(5,i)/U_av(1,i)*gamma/(gamma-1.0) + ( U_av(2,i)**2+U_av(3,i)**2 )/2.0 ) 
            
            !H= U(5)/U(1)*gamma / (gamma-1.0) + ( U(2)**2 + U(3)**2 ) /2.0 
           
        case default
      
            do j=1,5
                UL(j) = U(j,ncl) + psi(j,ncl) * dot_product( Grad(:,j,ncl),rL(:,i) )
                UR(j) = U(j,ncr) + psi(j,ncr) * dot_product( Grad(:,j,ncr),rR(:,i) )
                
               
            end do 
         
                
            nx=vector(1,i)/ds(i)
            ny=vector(2,i)/ds(i) 
    
            HL= UL(5)/UL(1)*gamma / (gamma-1.0) + ( UL(2)**2 + UL(3)**2 ) /2.0 
            HR= UR(5)/UR(1)*gamma / (gamma-1.0) + ( UR(2)**2 + UR(3)**2 ) /2.0 
    
            !seven Roe'S_SA averag_SAing_SA flow properities 
            rou_Roe= sqrt( UL(1)*UR(1) )
            u_Roe = ( UL(2)*sqrt(UL(1)) + UR(2)* sqrt(UR(1)) ) / ( sqrt(UL(1)) + sqrt(UR(1)) )
            v_Roe = ( UL(3)*sqrt(UL(1)) + UR(3)* sqrt(UR(1)) ) / ( sqrt(UL(1)) + sqrt(UR(1)) )
            H_Roe=  ( HL* sqrt(UL(1))  + HR*sqrt(UR(1)) ) / ( sqrt(UL(1)) + sqrt(UR(1)) )
    
            q2_Roe = u_Roe**2 + v_Roe**2
            c_Roe = sqrt( (gamma-1.0)*(H_Roe-q2_Roe/2.0) )
         
            Vn_Roe = u_Roe * nx + v_Roe * ny 
            !-------------------------------------------------
                
            VnL = UL(2)*nx + UL(3)*ny 
            VnR = UR(2)*nx + UR(3)*ny
    
            FcL(1) = UL(1)*VnL
            FCL(2) = UL(1)*UL(2)*VnL + nx*UL(5)
            FcL(3) = UL(1)*UL(3)*VnL + ny*UL(5)
            FcL(5) = UL(1)*HL*VnL
        
            FcR(1) = UR(1)*VnR   
            FCR(2) = UR(1)*UR(2)*VnR + nx*UR(5)
            FcR(3) = UR(1)*UR(3)*VnR + ny*UR(5)
            FcR(5) = UR(1)*HR*VnR
            
            
            !-------------------------------------------------------
            !HarT_SAen'S_SA entropy correction
            delta = 0.1*c_Roe  !o.1 times of the local sound speed
        
            lamda_c1 = abs(Vn_Roe - c_Roe)
            lamda_c234 = abs(Vn_Roe)
            lamda_c5 = abs(Vn_Roe + c_Roe)
        
            if( lamda_c1 .LE. delta )   lamda_c1 = (lamda_c1**2 + delta**2) / (2*delta)
            if( lamda_c234 .LE. delta ) lamda_c234 = (lamda_c234**2 + delta**2) / (2*delta)
            if( lamda_c5 .LE. delta )   lamda_c5 = (lamda_c5**2 + delta**2) / (2*delta)
            
    
            !-------------------------------------------------
         
            dF1(1) = lamda_c1*( UR(5) - UL(5) - rou_Roe*c_Roe*(VnR-VnL) ) / ( 2.0*c_Roe**2 )
            dF1(2) = lamda_c1*( UR(5) - UL(5) - rou_Roe*c_Roe*(VnR-VnL) ) / ( 2.0*c_Roe**2 ) * (u_Roe - c_Roe*nx)
            dF1(3) = lamda_c1*( UR(5) - UL(5) - rou_Roe*c_Roe*(VnR-VnL) ) / ( 2.0*c_Roe**2 ) * (v_Roe - c_Roe*ny)
            dF1(5) = lamda_c1*( UR(5) - UL(5) - rou_Roe*c_Roe*(VnR-VnL) ) / ( 2.0*c_Roe**2 ) * (H_Roe - c_Roe*Vn_Roe)
        
            dF5(1) = lamda_c5*( UR(5) - UL(5) + rou_Roe*c_Roe*(VnR-VnL) ) / ( 2.0*c_Roe**2 )  
            dF5(2) = lamda_c5*( UR(5) - UL(5) + rou_Roe*c_Roe*(VnR-VnL) ) / ( 2.0*c_Roe**2 ) * (u_Roe + c_Roe*nx)
            dF5(3) = lamda_c5*( UR(5) - UL(5) + rou_Roe*c_Roe*(VnR-VnL) ) / ( 2.0*c_Roe**2 ) * (v_Roe + c_Roe*ny)
            dF5(5) = lamda_c5*( UR(5) - UL(5) + rou_Roe*c_Roe*(VnR-VnL) ) / ( 2.0*c_Roe**2 ) * (H_Roe + c_Roe*Vn_Roe)
        
            dF234(1) = lamda_c234 * ( UR(1)-UL(1) - ( UR(5)-UL(5) ) / c_Roe**2 )
            dF234(2) = lamda_c234 * ( (UR(1)-UL(1) -( UR(5)-UL(5) ) / c_Roe**2) * u_Roe + rou_Roe*( UR(2)-UL(2) -(VnR-VnL)*nx ) )
            dF234(3) = lamda_c234 * ( (UR(1)-UL(1) -( UR(5)-UL(5) ) / c_Roe**2) * v_Roe + rou_Roe*( UR(3)-UL(3) -(VnR-VnL)*ny ) )
            dF234(5) = lamda_c234 * ( (UR(1)-UL(1) -( UR(5)-UL(5) ) / c_Roe**2) * q2_Roe/2.0 + rou_Roe*( u_Roe*(UR(2)-UL(2)) + v_Roe*(UR(3)-UL(3)) - Vn_Roe*(VnR-VnL) ) )
       
     
            Fc(:,ncl) = Fc(:,ncl) + 1.0/2 * ( FcL + FcR - dF1-dF234-dF5 )  * ds(i)
            Fc(:,ncr) = Fc(:,ncr) - 1.0/2 * ( FcL + FcR - dF1-dF234-dF5 )  * ds(i)
        
        end select
      
    end do
    
end subroutine


    


 
