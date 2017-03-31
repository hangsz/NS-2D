subroutine outputFreeFlow
    !-----------------------------------
        ! purpose: output the free inflow information
    !-----------------------------------
    implicit none
    real(8)::AOA(phase)
    integer::idFree
    integer::i
    !------------------------------------
    write(*,*)  "outputFreeFlow"
    do i=1,phase
        AoA(i) = att + att_ampl*sin(pi/2.0/50*itermax*i)
    end do
    
    open( newunit=idFree,file = "freeFlow/freeFlowInfo.dat")
        write(idFree,*)   att
        write(idFree,*)   phase
        write(idFree,*)   AOA
        write(idFree,*)  " ma  rou  p "
        write(idFree,*)   Ma_inf, rou_inf, P_inf
        
    close(idFree)
    
end subroutine 