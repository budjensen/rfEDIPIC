module colors
   use ParallelOperationValues
   use mpi
   implicit none

   ! Define colors for printing to the terminal
   character(len=5), parameter :: red = achar(27)//'[31m'        !> Red color
   character(len=5), parameter :: gre = achar(27)//'[32m'        !> Green color

   ! Define bold colors for printing to the terminal
   character(len=7), parameter :: bred = achar(27)//'[1;31m'     !> Bold red color
   character(len=7), parameter :: bgre = achar(27)//'[1;32m'     !> Bold green color

   character(len=4), parameter :: bld = achar(27)//'[1m'         !> Bold style

   ! Define the end of style command
   character(len=4), parameter :: normal = achar(27)//'[0m'      !> Normal style

end module colors

module testing
   use mpi
   use colors
   implicit none

contains
   !> Prints a message to the terminal if the test passed
   subroutine test_passed()
      print '(10x,A,"PASSED",A)', bgre, normal
   end subroutine test_passed

   !> Prints a message to the terminal if the test failed
   !! ierr is the number of errors reported by the test, if value is non-zero
   subroutine test_failed(ierr)
      implicit none
      integer, intent(in) :: ierr !> Error code

      if (ierr /= 0) then
         print '(10x,A,"FAILED",A," with",i4, " errors")', bred, normal, ierr
      else
         print '(10x,A,"FAILED",A)', bred, normal
      end if
   end subroutine test_failed

   !> Prints the test name and a test initiation message to the terminal.
   !! Optionally prints the rank of the process and the number of processes
   !! if the verbose flag is set to .true.
   subroutine test_header(name, rank, n_proc, verbose)
      integer, intent(in)          :: n_proc          !> Nmber of parallel processes
      integer, intent(in)          :: rank            !> Parallel rank
      character(len=*), intent(in) :: name            !> Name of test
      logical, intent(in)          :: verbose         !> Print verbose output
      integer                      :: ierr            !> Error code
   
      if (rank == 0) then
         print '(/2x"Running test:",3x,A,/)', trim(bld//name//normal)
      end if
      call mpi_barrier(MPI_COMM_WORLD, ierr)
   
      if (verbose .and. (rank > 0)) print '("  - running on client ",i1," of ",i1)', rank, n_proc - 1
      call mpi_barrier(MPI_COMM_WORLD, ierr)
   
   end subroutine test_header

   !> Prints the test name and a test completion message to the terminal.
   !! Optionally prints the rank of the process and the number of processes
   !! if the verbose flag is set to .true.
   subroutine test_footer(name, rank, n_proc, verbose)
      integer, intent(in)          :: n_proc          !> Nmber of parallel processes
      integer, intent(in)          :: rank            !> Parallel rank
      character(len=*), intent(in) :: name            !> Name of test
      logical, intent(in)          :: verbose         !> Print verbose output
      integer                      :: ierr            !> Error code
   
      if (rank == 0) then
         print '(/2x"Finished test:",3x,A,/)', trim(bld//name//normal)
      end if
      call mpi_barrier(MPI_COMM_WORLD, ierr)
   
      if (verbose .and. (rank > 0)) print '("  - on client ",i1," of ",i1)', rank, n_proc - 1
      call mpi_barrier(MPI_COMM_WORLD, ierr)
   
   end subroutine test_footer
end module testing