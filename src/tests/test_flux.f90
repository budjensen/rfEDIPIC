program test_flux
   use CurrentProblemValues, only : N_spec, N_part, N_cells, X_of_spec, VX_of_spec, Length, V_scl_ms, N_scl_m3
   use Diagnostics, only : flux_m2s1, N_of_probes
   use ParallelOperationValues, only : Rank_of_process, N_of_processes
   use colors
   use testing
   use mpi
   implicit none

   integer :: ierr, alloc_err, dealloc_err               ! Error codes
   integer :: node_idx, species_idx
   integer :: test_loc
   real(8) :: flux
   real(8) :: flux_at_node
   real(8) :: ans = 0.0                                  !! Expected result
   character(len=30) :: program_name = "test_flux"
   logical :: verbose = .false.

   ! Initialize MPI
   call mpi_init(ierr)
   call mpi_comm_rank(MPI_COMM_WORLD, Rank_of_process, ierr)
   call mpi_comm_size(MPI_COMM_WORLD, N_of_processes, ierr)

   if (N_of_processes == 1) then
      print '(/2x,"N_of_processes = 0. No processes available to run tests.")'
      stop
   end if

   ! Begin the test
   call test_header(program_name, Rank_of_process, N_of_processes, verbose)

   ! Set the common problem values
   V_scl_ms           = 1.0
   N_scl_m3           = 1.0
   N_cells            = 100
   N_part             = 1000
   N_of_probes        = 1
   Length             = N_part


   ! Allocate the arrays
   call CONFIGURE_PARTICLE_DYNAM_ARRAYS
   call AllocateDiagnosticArrays

   if (Rank_of_process /= 0) then
      ! Set the particle values for the clients
      X_of_spec(1)%part  = 5.0
      VX_of_spec(1)%part = 1.0
      X_of_spec(2)%part  = 5.0
      VX_of_spec(2)%part = -1.0
   else
      ! Set the particle values for the server
      X_of_spec(1)%part  = 5.0
      VX_of_spec(1)%part = 1.0
      X_of_spec(2)%part  = 5.0
      VX_of_spec(2)%part = -1.0
   end if

   ! Set the node index to 5
   node_idx = 5

   ! Verify each species passes
   do species_idx = 1, N_spec
      ! Call the flux_at_node function
      flux = flux_at_node(node_idx, species_idx)

      ! Calculate the expected result
      test_loc = int(X_of_spec(species_idx)%part(1))
      if (node_idx == test_loc) then
         ans = VX_of_spec(species_idx)%part(1) * dble(N_part(species_idx)) * dble(N_of_processes - 1)
      else
         ans = 0.0
      end if

      ! Print status
      if (Rank_of_process == 0) then
         print '(4x,"Verifying flux of species ",i1,"...")', species_idx
      end if

      ! Print the result
      if (Rank_of_process == 0) then
         if (flux /= ans) then
            call test_failed(0)
         else
            call test_passed()
         end if
      end if
   end do

   ! Now set node_idx = 4 and make sure the flux is zero
   node_idx = 4

   ! Verify each species passes
   do species_idx = 1, N_spec
      ! Call the flux_at_node function
      flux = flux_at_node(node_idx, species_idx)

      ! Calculate the expected result
      test_loc = int(X_of_spec(species_idx)%part(1))
      if (node_idx == test_loc) then
         ans = VX_of_spec(species_idx)%part(1) * dble(N_part(species_idx)) * dble(N_of_processes - 1)
      else
         ans = 0.0
      end if

      ! Print status
      if (Rank_of_process == 0) then
         print '(4x,"Verifying flux of species ",i1," is zero...")', species_idx
      end if

      ! Print the result
      if (Rank_of_process == 0) then
         if (flux /= ans) then
            call test_failed(0)
         else
            call test_passed()
         end if
      end if
   end do

   ! Now test the subroutine get_instantaneous_flux()
   ! Print status
   if (Rank_of_process == 0) then
      print '(4x,"Verifying flux across all nodes...")'
   end if

   ! Calculate the flux across all nodes
   call get_instantaneous_flux()

   ! Verify the results by making sure all elements are zero, except for node=5
   ierr = 0
   if (Rank_of_process > 0) then
      do species_idx = 1, N_spec
         do node_idx = 0, N_cells
            test_loc = int(X_of_spec(species_idx)%part(1))
            if (node_idx == test_loc) then
               ans = VX_of_spec(species_idx)%part(1) * dble(N_part(species_idx))
            else
               ans = 0.0
            end if
            if (flux_m2s1(node_idx, species_idx) /= ans) then
               print *, "flux_m2s1(", node_idx, ",", species_idx, ") = ", flux_m2s1(node_idx, species_idx), " /= ", ans
               ierr = ierr + 1
            end if
         end do
      end do
   end if

   ! Reduce the ierr to the master process
   call mpi_reduce(ierr, test_loc, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
   call mpi_barrier(MPI_COMM_WORLD, ierr)

   ! Print the result
   if (Rank_of_process == 0) then
      ierr = test_loc
      if (ierr > 0) then
         call test_failed(ierr)
      else
         call test_passed()
      end if
   end if

   ! Deallocate the arrays
   call REMOVE_PARTICLE_DYNAM_ARRAYS
   call DeallocateDiagnosticArrays

   ! Print the final message
   call test_footer(program_name, Rank_of_process, N_of_processes, verbose)

   ! Finalize MPI
   call mpi_finalize(ierr)

end program test_flux