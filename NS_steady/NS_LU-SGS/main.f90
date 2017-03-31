
program main
    use FVM
    implicit none

    !do
    !    write(*,*) "-------------------------------------------"
    !    write(*,*)  "is this simulation based on the new grid?"
    !    write(*,*)  "y-yes  n-no"
    !    read(*,*)   isNewGrid
    !    write(*,*) "-------------------------------------------"
    !    if( isNewGrid=="n" .OR. isNewGrid=="y")  exit
    !end do
    
    if( schemeNum=="1" .OR. schemeNum=="2")  then
        call Solver_SA
    else
        call Solver_SST
    end if
    
      
end program
