subroutine random_init
USE CurrentProblemValues, ONLY : seed
integer :: values(1:8), k, j
!integer, dimension(:), allocatable :: seed

call date_and_time(values=values)
call random_seed(size=k)
allocate(seed(1:k))
j = values(8)**5 - 1
if (j.le.0) j = (2**31 - 1) + j !** if integer overflow occurs
seed(:) = j
call random_seed(put=seed)
write(*,*) 'Initializing pseudorandom numbers:'
write(*,*)   ((seed(j)), j = 1,k)
write(*,*) ((values(j)), j = 1,8)

return
end subroutine random_init
  
