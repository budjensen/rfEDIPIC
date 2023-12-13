!===========================================
program MainProg

   use mpi
   use CurrentProblemValues
   use ParallelOperationValues

   implicit none

!  include 'mpif.h'

   integer ierr
   real(8) start, finish

   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, Rank_of_process, ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, N_of_processes, ierr)

!  print *, 'my rank is ', Rank_of_process

   call INITIATE_PARAMETERS         ! all, with some differences, distribute particles over nodes here

   call INITIATE_MC_COLLISIONS
   call INITIATE_SE_EMISSION        ! all, with some differences
   call INITIATE_ELECTRON_INJECTION
   call INITIATE_BEAM_IN_PLASMA
   call INITIATE_COULOMB_COLLISIONS

   call INITIATE_DIAGNOSTICS

   if (Rank_of_process.EQ.0) then
      call PREPARE_TIME_DEPENDENCE_DATAFILES         ! server only
!restore this later     call INITIATE_TEST_PARTICLES                   ! server only
   end if

   start = MPI_WTIME()

   do T_cntr = Start_T_cntr, Max_T_cntr

      if (Rank_of_process.eq.0) start = MPI_WTIME()

!!! print *, T_cntr, Rank_of_process, N_part(1), N_part(2)

      if (T_cntr.NE.Start_T_cntr) then
         call CALCULATE_STR_CHARGE_DENSITY                                   ! server and clients, with difference
         call CALCULATE_STR_LONG_ELECTR_FIELD                                ! server and clients, with difference

         if (Rank_of_process.GT.0) then
            call FINAL_PUSH                                                  ! clients
            call INJECT_ADDITIONAL_EI_PAIRS_LEFT_WALL                        ! clients
            call INJECT_ELECTRONS_AT_WALLS                                   ! clients
         else
!restore this later           call TEST_PARTICLES_COMBINED_PUSH_AND_SAVE                       ! server
         end if

         call MANAGE_WALL_CHARGE_DENSITIES                                   ! server and clients, with difference

         call PROCESS_COLLISIONS_WITH_ATOMS                                  ! server and clients, with difference

         call PROCESS_COLLISIONS_EL_EL_LNGV                                  ! server and clients, with difference

         call START_BEAM_IN_PLASMA                                           ! server and clients, with difference, only one time

         if (Rank_of_process.GT.0) call REARRANGE_ACTIVE_BLOCKS              ! clients

!        call CLEAR_TAGS                                 ! clients

         if (Rank_of_process.GT.0) call CLEAR_ELECTRON_LIST                  ! clients
         if (Rank_of_process.GT.0) call CLEAR_ION_LIST                       ! clients

         call MPI_BARRIER(MPI_COMM_WORLD, ierr)                              ! all, synchronization   !!!! comment me !!!!

      end if

      ! Save diagnostics data, if inside of a diagnostics window [Start_diag_Tcntr, Finish_diag_Tcntr]
      call do_DIAGNOSTICS_STEP_1                                             ! server and clients, with differences

      call CREATE_SNAPSHOT                                                   ! server and clients, with differences

      ! If we reach Finish_diag_Tcntr, clear diagnostics arrays and set new a new diagnostics window
      call do_DIAGNOSTICS_STEP_2                                             ! server and clients, with differences

      if (SaveCheck_step.GE.0) call SAVE_CHECKPOINT                          ! server and clients, with differences

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)                                 ! all, synchronization

      if (Rank_of_process.GT.0) then
         call PREALLOCATE_ELECTRON_LIST                                      ! clients
         call PREALLOCATE_ION_LIST                                           ! clients
         call PREDICTING_PUSH                                                ! clients
      else
!restore this later        call ACTIVATE_TEST_PARTICLE_RUN                                     ! server
      end if

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)                                 ! all, synchronization

!if (Rank_of_process.eq.0) then
!   finish = MPI_WTIME()
!   if (mod(T_cntr,17).eq.0) print '(2x,"**** cycle ",i8," time is  : ", f11.6," sec")', T_cntr, finish - start
!end if

   end do

   finish = MPI_WTIME()

   print '(2x,"**** Process ",i3" : Simulation time is  : ", f11.3," sec")', Rank_of_process, finish - start

   if (Rank_of_process.GT.0) call CLEAR_ELECTRON_LIST                        ! clients
   if (Rank_of_process.GT.0) call CLEAR_ION_LIST                             ! clients

   call REMOVE_MESH_ARRAYS                               ! all
   if (Rank_of_process.GT.0) call REMOVE_PARTICLE_DYNAM_ARRAYS               ! clients (server calls it in INITIATE_PARAMETERS : INITIATE_PARTICLE_DYNAM_ARRAYS)

   call REMOVE_CRSECT_ARRAYS
   call REMOVE_COLLSPEC_ARRAYS
   call REMOVE_COLLISION_ARRAYS
   call REMOVE_LANGEVIN_ARRAYS

   call FINISH_SNAPSHOTS
   call FINISH_DIAGNOSTICS
!restore this later  if (Rank_of_process.EQ.0) call REMOVE_TEST_PARTICLE_ARRAYS                ! server only

   call MPI_FINALIZE(ierr)

end program MainProg


