subroutine flowInit_dynamic
    implicit none
    
    integer::idRead
    character( len = 200)::filename
    logical::alive
    
    !----------------------------------------------------
    filename = trim(caseRoute)//"/"//trim(gridNum)//"/"//trim(fileNum)//"/steady/steady.dat"
    inquire(file =filename,exist =alive)
    
    if ( .NOT. alive) then
        write(*,*)  filename,'not exist.'
        stop   
    end if
    
    select case( turModelNum )
    case( 0 )
        call flowInit_turNone_dynamic(filename)
    case( 1 )
        call flowInit_SA_dynamic(filename)
    case( 2 )
        call flowInit_SST_dynamic(filename)    
    case default
        
    end select
    
end subroutine
 
    
subroutine flowInit_turNone_dynamic(filename)

    implicit none
    character( len=200),intent(in)::filename
    integer::i
    integer::idRead
    
    open(newunit=idRead,file=filename)
    
        read(idRead,*) 
        read(idRead,*) 
        read(idRead,*) 
        read(idRead,*) 
        read(idRead,*) 
     
        do i=1,nnodes
            read(idRead,*) !xy(1,i)
        end do
        do i=1,nnodes
            read(idRead,*) !xy(2,i)
        end do
    
        do i=1,ncells
            read(idRead,*) U(5,i)
        end do
        do i=1,ncells
            read(idRead,*) U(1,i)
        end do
        do i=1,ncells
            read(idRead,*) U(2,i)
        end do
        do i=1,ncells
            read(idRead,*) U(3,i)
        end do
    
        do i=1,ncells
            read(idRead,*) !icell(:,i)
        end do
    
    close(idRead)
    
    a_inf=sqrt(p_inf*gamma/rou_inf)
    u_inf=Ma_inf*a_inf*cosd(att)        
    v_inf=Ma_inf*a_inf*sind(att)
   
    W(1,:)=U(1,:)
    W(2,:)=U(1,:)*U(2,:)
    W(3,:)=U(1,:)*U(3,:)
    W(5,:)=U(5,:)/(gamma-1.0) + U(1,:)*( U(2,:)**2 + U(3,:)**2 )/2.0   
    
    muT =0.0
    muT_av = 0.0

end subroutine
    
    
subroutine flowInit_SA_dynamic(filename)      !initialize the flow field
   ! purpose: read the stady flow properties
  
    implicit none 
    character( len = 200),intent(in)::filename
    integer::i
 
    integer::idRead
   
    open(newunit=idRead,file=filename)
    
        read(idRead,*) 
        read(idRead,*) 
        read(idRead,*) 
        read(idRead,*) 
        read(idRead,*) 
     
        do i=1,nnodes
            read(idRead,*) !xy(1,i)
        end do
        do i=1,nnodes
            read(idRead,*) !xy(2,i)
        end do
    
        do i=1,ncells
            read(idRead,*) U(5,i)
        end do
        do i=1,ncells
            read(idRead,*) U(1,i)
        end do
        do i=1,ncells
            read(idRead,*) U(2,i)
        end do
        do i=1,ncells
            read(idRead,*) U(3,i)
        end do
        do i=1,ncells
            read(idRead,*) nu_SA(i)
        end do
        do i=1,ncells
            read(idRead,*) !icell(:,i)
        end do
    
    close(idRead)
    
    a_inf=sqrt(p_inf*gamma/rou_inf)
    u_inf=Ma_inf*a_inf*cosd(att)        
    v_inf=Ma_inf*a_inf*sind(att)
   
    W(1,:)=U(1,:)
    W(2,:)=U(1,:)*U(2,:)
    W(3,:)=U(1,:)*U(3,:)
    W(5,:)=U(5,:)/(gamma-1.0) + U(1,:)*( U(2,:)**2 + U(3,:)**2 )/2.0  
    
    !initilize the turbulent flow
    ! nuT = 0.1*muL_inf / rou_inf
    W_SA = U(1,:) * nu_SA
    
   
end subroutine

subroutine flowInit_SST_dynamic(filename)      !initialize the flow field

    implicit none
    character( len = 200),intent(in)::filename
    integer::i
    integer::idRead
    ! purpose: read the stady flow properties
  
    open(newunit=idRead,file=filename)
    
    read(idRead,*) 
    read(idRead,*) 
    read(idRead,*) 
    read(idRead,*) 
    read(idRead,*) 
     
    do i=1,nnodes
        read(idRead,*) !xy(1,i)
    end do
    do i=1,nnodes
        read(idRead,*) !xy(2,i)
    end do
    
    do i=1,ncells
        read(idRead,*) U(5,i)
    end do
    do i=1,ncells
        read(idRead,*) U(1,i)
    end do
    do i=1,ncells
        read(idRead,*) U(2,i)
    end do
    do i=1,ncells
        read(idRead,*) U(3,i)
    end do
    do i=1,ncells
        read(idRead,*) muT(i)
    end do
    do i=1,ncells
        read(idRead,*) U_SST(1,i)
    end do
    do i=1,ncells
        read(idRead,*) U_SST(2,i)
    end do
    do i=1,ncells
        read(idRead,*) !icell(:,i)
    end do
    
    close(idRead)
   
    a_inf=sqrt(p_inf*gamma/rou_inf)
    u_inf=Ma_inf*a_inf*cosd(att)        
    v_inf=Ma_inf*a_inf*sind(att)
    
    W(1,:)=U(1,:)
    W(2,:)=U(1,:)*U(2,:)
    W(3,:)=U(1,:)*U(3,:)
    W(5,:)=U(5,:)/(gamma-1.0) + U(1,:)*(U(2,:)**2 + U(3,:)**2)/2.0  + U(1,:)*U_SST(1,:)
    
    !initilize the turbulent flow
    W_SST(1,:) = U(1,:) * U_SST(1,:)
    W_SST(2,:) = U(1,:) * U_SST(2,:)
   
    
    
end subroutine
