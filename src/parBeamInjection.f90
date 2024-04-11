
SUBROUTINE INITIATE_INJECTION

   use mpi
   USE ParallelOperationValues
   USE ELectronInjection
   USE CurrentProblemValues, ONLY : T_e_eV, BC_flag, N_box_vel, N_max_vel, r_max_vel, N_of_particles_cell, delta_t_s, e_Cl, N_plasma_m3, delta_x_m, Qs, Ms
   IMPLICIT NONE

!  INCLUDE 'mpif.h'

   LOGICAL exists
   REAL(8) length_of_injection
   INTEGER left_for_injection

   CHARACTER (77) buf
   INTEGER ierr

   REAL(8) temp

   INQUIRE (FILE = 'ssc_einject.dat', EXIST = exists)
   CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

   OPEN (9, FILE = 'ssc_einject.dat')
   IF(exists) THEN

      IF (Rank_of_process.EQ.0) &
      & PRINT '(2x,"Process ",i3," : The electron injection data file is found. Reading the data file...")', Rank_of_process

      READ (9, '(A77)') buf ! ****************** LEFT WALL, CONSTANT ELECTRON INJECTION *******************")')
      READ (9, '(A77)') buf ! ---           0 = turned off                                              ---")')
      READ (9, '(A77)') buf ! ---           1 = monoenergetic beam                                      ---")')
      READ (9, '(A77)') buf ! -------d----- 2 = warm maxwellian source ------------------------------------")')
      READ (9, '(7x,i1)') eBeamInjectFlag_left
      READ (9, '(A77)') buf ! --dddddd.ddd- Start injection at this moment, [ns] (timesteps if < 0) -------")')
      READ (9, '(2x,f10.3)') Delay_of_e_injection_ns_left

      READ (9, '(A77)') buf ! ====================== ELECTRON INJECTION, PARAMETERS =======================")')
      READ (9, '(A77)') buf ! --d.ddddE#dd- Energy (cold) / Temperature (warm), [eV] ----------------------")')
      READ (9, '(1x,e11.4)') eBeam_energy_eV_left
      READ (9, '(A77)') buf ! --d.ddddE#dd- Injection current magnitude [A/m^2] ---------------------------")')
      READ (9, '(1x,e11.4)') eBeam_J_Am2_left

      READ (9, '(A77)') buf ! ********************* LEFT WALL, CONSTANT ION INJECTION *********************")')
      READ (9, '(A77)') buf ! ---           0 = turned off                                              ---")')
      READ (9, '(A77)') buf ! ---           1 = monoenergetic beam                                      ---")')
      READ (9, '(A77)') buf ! -------d----- 2 = warm maxwellian source ------------------------------------")')
      READ (9, '(7x,i1)') iBeamInjectFlag_left
      READ (9, '(A77)') buf ! --dddddd.ddd- Start injection at this moment, [ns] (timesteps if < 0) -------")')
      READ (9, '(2x,f10.3)') Delay_of_i_injection_ns_left

      READ (9, '(A77)') buf ! ========================= ION INJECTION, PARAMETERS =========================")')
      READ (9, '(A77)') buf ! --d.ddddE#dd- Energy (cold) / Temperature (warm), [eV] ----------------------")')
      READ (9, '(1x,e11.4)') iBeam_energy_eV_left
      READ (9, '(A77)') buf ! --d.ddddE#dd- Injection current magnitude [A/m^2] ---------------------------")')
      READ (9, '(1x,e11.4)') iBeam_J_Am2_left

      READ (9, '(A77)') buf ! ****************** RIGHT WALL, CONSTANT ELECTRON INJECTION ******************")')
      READ (9, '(A77)') buf ! ---           0 = turned off                                              ---")')
      READ (9, '(A77)') buf ! ---           1 = monoenergetic beam                                      ---")')
      READ (9, '(A77)') buf ! -------d----- 2 = warm maxwellian source ------------------------------------")')
      READ (9, '(7x,i1)') eBeamInjectFlag_right
      READ (9, '(A77)') buf ! --dddddd.ddd- Start injection at this moment, [ns] (timesteps if < 0) -------")')
      READ (9, '(2x,f10.3)') Delay_of_e_injection_ns_right

      READ (9, '(A77)') buf ! ====================== ELECTRON INJECTION, PARAMETERS =======================")')
      READ (9, '(A77)') buf ! --d.ddddE#dd- Energy (cold) / Temperature (warm), [eV] ----------------------")')
      READ (9, '(1x,e11.4)') eBeam_energy_eV_right
      READ (9, '(A77)') buf ! --d.ddddE#dd- Injection current magnitude [A/m^2] ---------------------------")')
      READ (9, '(1x,e11.4)') eBeam_J_Am2_right

      READ (9, '(A77)') buf ! ******************** RIGHT WALL, CONSTANT ION INJECTION *********************")')
      READ (9, '(A77)') buf ! ---           0 = turned off                                              ---")')
      READ (9, '(A77)') buf ! ---           1 = monoenergetic beam                                      ---")')
      READ (9, '(A77)') buf ! -------d----- 2 = warm maxwellian source ------------------------------------")')
      READ (9, '(7x,i1)') iBeamInjectFlag_right
      READ (9, '(A77)') buf ! --dddddd.ddd- Start injection at this moment, [ns] (timesteps if < 0) -------")')
      READ (9, '(2x,f10.3)') Delay_of_i_injection_ns_right

      READ (9, '(A77)') buf ! ========================= ION INJECTION, PARAMETERS =========================")')
      READ (9, '(A77)') buf ! --d.ddddE#dd- Energy (cold) / Temperature (warm), [eV] ----------------------")')
      READ (9, '(1x,e11.4)') iBeam_energy_eV_right
      READ (9, '(A77)') buf ! --d.ddddE#dd- Injection current magnitude [A/m^2] ---------------------------")')
      READ (9, '(1x,e11.4)') iBeam_J_Am2_right

      READ (9, '(A77)') buf ! ******************************* SMART TAGS **********************************")')
      READ (9, '(A77)') buf ! ---           0 = turned off (changing velocity sign does not change tag) ---")')
      READ (9, '(A77)') buf ! -------d----- 1 = turned on -------------------------------------------------")')
      READ (9, '(7x,i1)') UseSmartTagsFlag

   ELSE

      eBeamInjectFlag_left          = 0           ! injection at the left wall is turned off
      Delay_of_e_injection_ns_left  = 0.0_8
      eBeam_energy_eV_left          = 0.0_8
      eBeam_J_Am2_left              = 0.0_8

      iBeamInjectFlag_left          = 0           ! injection at the left wall is turned off
      Delay_of_i_injection_ns_left  = 0.0_8
      iBeam_energy_eV_left          = 0.0_8
      iBeam_J_Am2_left              = 0.0_8

      eBeamInjectFlag_right         = 0           ! injection at the right wall is turned off
      Delay_of_e_injection_ns_right = 0.0_8
      eBeam_energy_eV_right         = 0.0_8
      eBeam_J_Am2_right             = 0.0_8

      iBeamInjectFlag_right         = 0           ! injection at the right wall is turned off
      Delay_of_i_injection_ns_right = 0.0_8
      iBeam_energy_eV_right         = 0.0_8
      iBeam_J_Am2_right             = 0.0_8

      UseSmartTagsFlag = 1                      ! smart tags are on

      IF (Rank_of_process.EQ.0) THEN

         PRINT '(/2x,"Process ",i3," : File with the name ssc_einject.dat not found. Use the default settings ...")'

         PRINT '(2x,"Process ",i3," : Create ssc_beaminject.dat file . . .")', Rank_of_process

         WRITE (9, '("****************** LEFT WALL, CONSTANT ELECTRON INJECTION *******************")')
         WRITE (9, '("---           0 = turned off                                              ---")')
         WRITE (9, '("---           1 = monoenergetic beam                                      ---")')
         WRITE (9, '("-------d----- 2 = warm maxwellian source ------------------------------------")')
         WRITE (9, '(7x,i1)') eBeamInjectFlag_left
         WRITE (9, '("--dddddd.ddd- Start injection at this moment, [ns] (timesteps if < 0) -------")')
         WRITE (9, '(2x,f10.3)') Delay_of_e_injection_ns_left

         WRITE (9, '("====================== ELECTRON INJECTION, PARAMETERS =======================")')
         WRITE (9, '("--d.ddddE#dd- Energy (cold) / Temperature (warm), [eV] ----------------------")')
         WRITE (9, '(1x,e11.4)') eBeam_energy_eV_left
         WRITE (9, '("--d.ddddE#dd- Injection current magnitude [A/m^2] ---------------------------")')
         WRITE (9, '(1x,e11.4)') eBeam_J_Am2_left

         WRITE (9, '("********************* LEFT WALL, CONSTANT ION INJECTION *********************")')
         WRITE (9, '("---           0 = turned off                                              ---")')
         WRITE (9, '("---           1 = monoenergetic beam                                      ---")')
         WRITE (9, '("-------d----- 2 = warm maxwellian source ------------------------------------")')
         WRITE (9, '(7x,i1)') iBeamInjectFlag_left
         WRITE (9, '("--dddddd.ddd- Start injection at this moment, [ns] (timesteps if < 0) -------")')
         WRITE (9, '(2x,f10.3)') Delay_of_i_injection_ns_left
   
         WRITE (9, '("========================= ION INJECTION, PARAMETERS =========================")')
         WRITE (9, '("--d.ddddE#dd- Energy (cold) / Temperature (warm), [eV] ----------------------")')
         WRITE (9, '(1x,e11.4)') iBeam_energy_eV_left
         WRITE (9, '("--d.ddddE#dd- Injection current magnitude [A/m^2] ---------------------------")')
         WRITE (9, '(1x,e11.4)') iBeam_J_Am2_left

         WRITE (9, '("****************** RIGHT WALL, CONSTANT ELECTRON INJECTION ******************")')
         WRITE (9, '("---           0 = turned off                                              ---")')
         WRITE (9, '("---           1 = monoenergetic beam                                      ---")')
         WRITE (9, '("-------d----- 2 = warm maxwellian source ------------------------------------")')
         WRITE (9, '(7x,i1)') eBeamInjectFlag_right
         WRITE (9, '("--dddddd.ddd- Start injection at this moment, [ns] (timesteps if < 0) -------")')
         WRITE (9, '(2x,f10.3)') Delay_of_e_injection_ns_right

         WRITE (9, '("====================== ELECTRON INJECTION, PARAMETERS =======================")')
         WRITE (9, '("--d.ddddE#dd- Energy (cold) / Temperature (warm), [eV] ----------------------")')
         WRITE (9, '(1x,e11.4)') eBeam_energy_eV_right
         WRITE (9, '("--d.ddddE#dd- Injection current magnitude [A/m^2] ---------------------------")')
         WRITE (9, '(1x,e11.4)') eBeam_J_Am2_right

         WRITE (9, '("********************* RIGHT WALL, CONSTANT ION INJECTION ********************")')
         WRITE (9, '("---           0 = turned off                                              ---")')
         WRITE (9, '("---           1 = monoenergetic beam                                      ---")')
         WRITE (9, '("-------d----- 2 = warm maxwellian source ------------------------------------")')
         WRITE (9, '(7x,i1)') iBeamInjectFlag_right
         WRITE (9, '("--dddddd.ddd- Start injection at this moment, [ns] (timesteps if < 0) -------")')
         WRITE (9, '(2x,f10.3)') Delay_of_i_injection_ns_right
   
         WRITE (9, '("========================= ION INJECTION, PARAMETERS =========================")')
         WRITE (9, '("--d.ddddE#dd- Energy (cold) / Temperature (warm), [eV] ----------------------")')
         WRITE (9, '(1x,e11.4)') iBeam_energy_eV_right
         WRITE (9, '("--d.ddddE#dd- Injection current magnitude [A/m^2] ---------------------------")')
         WRITE (9, '(1x,e11.4)') iBeam_J_Am2_right

         WRITE (9, '("******************************* SMART TAGS **********************************")')
         WRITE (9, '("---           0 = turned off (changing velocity sign does not change tag) ---")')
         WRITE (9, '("-------d----- 1 = turned on -------------------------------------------------")')
         WRITE (9, '(7x,i1)') UseSmartTagsFlag

      END IF

   END IF

   CLOSE (9, STATUS = 'KEEP')

! for the left wall electrons ----------------------------
   IF (eBeamInjectFlag_left.EQ.0) THEN
      IF (Rank_of_process.EQ.0) PRINT '(/2x,"The electron injection at the left wall is turned off ...")'
      inject_e_every_this_many_timesteps_left = 1

   ELSE
      temp = (eBeam_J_Am2_left / e_Cl) / (N_plasma_m3 * delta_x_m / (delta_t_s * DBLE(N_of_particles_cell)))
      N_e_to_inject_total_left = INT(temp)

      IF (N_e_to_inject_total_left.GT.0) THEN
         const_N_e_to_inject_by_proc_left = N_e_to_inject_total_left / (N_of_processes-1)     ! each node will emit at least this number of electrons
         variable_N_e_to_inject_left = MOD(N_e_to_inject_total_left, N_of_processes-1)        ! leftover, will be distributed between all client nodes
         inject_e_every_this_many_timesteps_left = 1

      ELSE IF (temp.GE.0.001_8) THEN  ! let the minimal permitted current is when a macroparticle is injected every 1000 timesteps
         const_N_e_to_inject_by_proc_left = 0                               ! only one macroparticle will be injected by one (each time different) process
         variable_N_e_to_inject_left = 1                                    ! every inject_every_this_many_timesteps_left timesteps
         inject_e_every_this_many_timesteps_left = INT(1.0_8 / temp) + 1    ! note that skip_timesteps_beam_inject_left >= 2

      ELSE
         const_N_e_to_inject_by_proc_left = 0          ! if the beam was requested but its current is too small
         variable_N_e_to_inject_left = 0               ! nothing will be injected
         inject_e_every_this_many_timesteps_left = 1   !

      END IF

      IF (Rank_of_process.GT.0) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
         PRINT '(2x,"Process ",i3," will inject at least ",i5," electron macroparticles at the LEFT boundary EACH timestep")', &
         & Rank_of_process, const_N_e_to_inject_by_proc_left
      ELSE
         ! for the server process print out the total number of macroparticles to be injected
         PRINT '(/2x,"Constant REQUESTED flux of electron injection from the LEFT WALL is ",e14.7," m^-2 s^-1")', (eBeam_J_Am2_left / e_Cl)
         PRINT '(2x,"TOTAL INTEGER number of electron macroparticles to be injected from the LEFT boundary at each timestep is : ",i7)', N_e_to_inject_total_left
         PRINT '(2x,"TOTAL REAL    number of electron macroparticles to be injected from the LEFT boundary at each timestep is : ",e12.5)', temp
         PRINT '(2x,"COMMON injection leftover distributed between all nodes is : ",i5)', variable_N_e_to_inject_left
         PRINT '(2x,"The leftover macroparticle(s) will be injected every ",i5," timestep(s)")', inject_e_every_this_many_timesteps_left
         PRINT '("-------")'
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)   ! place the synchronization point here to print server's output before the clients'
      END IF

      VX_e_beam_left = SQRT(eBeam_energy_eV_left / T_e_eV) * N_box_vel

      ! length_of_injection = SQRT(Beam_energy_eV_left / T_e_eV) / N_max_vel
      length_of_injection = SQRT(eBeam_energy_eV_left / T_e_eV) / r_max_vel
      IF (N_e_to_inject_total_left.GT.0) THEN
         Gap_betw_e_part_left = length_of_injection / N_e_to_inject_total_left
      ELSE
         Gap_betw_e_part_left = 0.0_8
      END IF

      ! if the delay of the start of injection is negative then it is given in timesteps, assumed to be integer !!!
      IF (Delay_of_e_injection_ns_left.LT.0.0_8) THEN
         IF (Rank_of_process.EQ.0) PRINT &
         &'(/2x,"The delay of electron injection is given in timesteps (negative value of Delay_of_e_injection_ns_left was read from the file)")'
         Delay_of_e_injection_T_cntr_left = - Delay_of_e_injection_ns_left
         Delay_of_e_injection_ns_left = Delay_of_e_injection_T_cntr_left * delta_t_s * 1.0e9

         IF (Rank_of_process.EQ.0) PRINT '(/2x,"The recalculated delay of electron injection is : ", f10.3)', Delay_of_e_injection_ns_left

      ELSE
         Delay_of_e_injection_T_cntr_left  = 1.0d-9 * Delay_of_e_injection_ns_left / delta_t_s                               ! Number of time steps
      END IF

   END IF

! for the left wall ions ----------------------------
   if (iBeamInjectFlag_left.eq.0) then
      if (Rank_of_process.eq.0) print '(/2x,"The ion injection at the left wall is turned off ...")'
      inject_i_every_this_many_timesteps_left = 1

   else
      temp = (iBeam_J_Am2_left / (Qs(2) * e_Cl)) / (N_plasma_m3 * delta_x_m / (delta_t_s * dble(N_of_particles_cell)))
      N_i_to_inject_total_left = int(temp)

      if (N_i_to_inject_total_left.gt.0) then
         const_N_i_to_inject_by_proc_left = N_i_to_inject_total_left / (N_of_processes-1)     ! each node will emit at least this number of ions
         variable_N_i_to_inject_left = mod(N_i_to_inject_total_left, N_of_processes-1)        ! leftover, will be distributed between all client nodes
         inject_i_every_this_many_timesteps_left = 1

      else if (temp.ge.0.001_8) then  ! minimal permitted current is one macroparticle injected every 1000 timesteps
         const_N_i_to_inject_by_proc_left = 0                               ! only one macroparticle will be injected by one (each time different) process
         variable_N_i_to_inject_left = 1                                    ! every inject_every_this_many_timesteps_left timesteps
         inject_i_every_this_many_timesteps_left = int(1.0_8 / temp) + 1    ! note that skip_timesteps_beam_inject_left >= 2

      else
         const_N_i_to_inject_by_proc_left = 0          ! if the beam was requested but its current is too small
         variable_N_i_to_inject_left = 0               ! nothing will be injected
         inject_i_every_this_many_timesteps_left = 1   !

      end if

      if (Rank_of_process.gt.0) then
         call mpi_barrier(MPI_COMM_WORLD, ierr)
         print '(2x,"Process ",i3," will inject at least ",i5," ion macroparticles at the LEFT boundary EACH timestep")', &
         & Rank_of_process, const_N_e_to_inject_by_proc_left
      else
         ! for the server process print out the total number of macroparticles to be injected
         print '(/2x,"Constant REQUESTED flux of ion injection from the LEFT WALL is ",e14.7," m^-2 s^-1")', (iBeam_J_Am2_left / e_Cl)
         print '(2x,"TOTAL INTEGER number of ion macroparticles to be injected from the LEFT boundary at each timestep is : ",i7)', N_i_to_inject_total_left
         print '(2x,"TOTAL REAL    number of ion macroparticles to be injected from the LEFT boundary at each timestep is : ",e12.5)', temp
         print '(2x,"COMMON injection leftover distributed between all nodes is : ",i5)', variable_N_i_to_inject_left
         print '(2x,"The leftover macroparticle(s) will be injected every ",i5," timestep(s)")', inject_i_every_this_many_timesteps_left
         print '("-------")'
         call mpi_barrier(MPI_COMM_WORLD, ierr)   ! place the synchronization point here to print server's output before the clients'
      end if

      VX_i_beam_left = sqrt(iBeam_energy_eV_left / (T_e_eV * Ms(2))) * N_box_vel

      ! if the delay of the start of injection is negative then it is given in timesteps, assumed to be integer !!!
      if (Delay_of_i_injection_ns_left.lt.0.0_8) then
         if (Rank_of_process.eq.0) print &
         &'(/2x,"The delay of ion injection is given in timesteps (negative value of Delay_of_i_injection_ns_left was read from the file)")'
         Delay_of_i_injection_T_cntr_left = - Delay_of_i_injection_ns_left
         Delay_of_i_injection_ns_left = Delay_of_i_injection_T_cntr_left * delta_t_s * 1.0e9

         if (Rank_of_process.eq.0) print '(/2x,"The recalculated delay of ion injection is : ", f10.3)', Delay_of_i_injection_ns_left

      else
         Delay_of_i_injection_T_cntr_left  = 1.0d-9 * Delay_of_i_injection_ns_left / delta_t_s                               ! Number of time steps
      end if

   end if

! for the right wall electrons -------------------------
   IF (eBeamInjectFlag_right.EQ.0) THEN
      IF (Rank_of_process.EQ.0) PRINT '(/2x,"The electron injection at the right wall is turned off ...")'
      inject_e_every_this_many_timesteps_right = 1

   ELSE
      temp = (eBeam_J_Am2_right / e_Cl) / (N_plasma_m3 * delta_x_m / (delta_t_s * DBLE(N_of_particles_cell)))
      N_e_to_inject_total_right = INT(temp)

      IF (N_e_to_inject_total_right.GT.0) THEN
         const_N_e_to_inject_by_proc_right = N_e_to_inject_total_right / (N_of_processes-1)     ! each node will emit at least this number of electrons
         variable_N_e_to_inject_right = MOD(N_e_to_inject_total_right, N_of_processes-1)        ! leftover, will be distributed between all client nodes
         inject_e_every_this_many_timesteps_right = 1

      ELSE IF (temp.GE.0.001_8) THEN  ! let the minimal permitted current is when a macroparticle is injected every 1000 timesteps
         const_N_e_to_inject_by_proc_right = 0                               ! only one macroparticle will be injected by one (each time different) process
         variable_N_e_to_inject_right = 1                                    ! every inject_every_this_many_timesteps_right timesteps
         inject_e_every_this_many_timesteps_right = INT(1.0_8 / temp) + 1    ! note that skip_timesteps_beam_inject_right >= 2

      ELSE
         const_N_e_to_inject_by_proc_right = 0          ! if the beam was requested but its current is too small
         variable_N_e_to_inject_right = 0               ! nothing will be injected
         inject_e_every_this_many_timesteps_right = 1   !

      END IF

      IF (Rank_of_process.GT.0) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
         PRINT '(2x,"Process ",i3," will inject at least ",i5," electron macroparticles at the RIGHT boundary EACH timestep")', &
         & Rank_of_process, const_N_e_to_inject_by_proc_right
      ELSE
         ! for the server process print out the total number of macroparticles to be injected
         PRINT '(/2x,"Constant REQUESTED flux of electron injection from the RIGHT WALL is ",e14.7," m^-2 s^-1")', (eBeam_J_Am2_right / e_Cl)
         PRINT '(2x,"TOTAL INTEGER number of electron macroparticles to be injected from the RIGHT boundary at each timestep is : ",i7)', N_e_to_inject_total_right
         PRINT '(2x,"TOTAL REAL    number of electron macroparticles to be injected from the RIGHT boundary at each timestep is : ",e12.5)', temp
         PRINT '(2x,"COMMON injection leftover distributed between all nodes is : ",i5)', variable_N_e_to_inject_right
         PRINT '(2x,"The leftover macroparticle(s) will be injected every ",i5," timestep(s)")', inject_e_every_this_many_timesteps_right
         PRINT '("-------")'
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)   ! place the synchronization point here to print server's output before the clients'
      END IF

      VX_e_beam_right = SQRT(eBeam_energy_eV_right / T_e_eV) * N_box_vel

      ! length_of_injection = SQRT(Beam_energy_eV_right / T_e_eV) / N_max_vel
      length_of_injection = SQRT(eBeam_energy_eV_right / T_e_eV) / r_max_vel
      IF (N_e_to_inject_total_right.GT.0) THEN
         Gap_betw_e_part_right = length_of_injection / N_e_to_inject_total_right
      ELSE
         Gap_betw_e_part_right = 0.0_8
      END IF

      ! if the delay of the start of injection is negative then it is given in timesteps, assumed to be integer !!!
      IF (Delay_of_e_injection_ns_right.LT.0.0_8) THEN
         IF (Rank_of_process.EQ.0) PRINT &
         &'(/2x,"The delay of electron injection is given in timesteps (negative value of Delay_of_e_injection_ns_right was read from the file)")'
         Delay_of_e_injection_T_cntr_right = - Delay_of_e_injection_ns_right
         Delay_of_e_injection_ns_right = Delay_of_e_injection_T_cntr_right * delta_t_s * 1.0e9

         IF (Rank_of_process.EQ.0) PRINT '(/2x,"The recalculated delay of electron injection is : ", f10.3)', Delay_of_e_injection_ns_right

      ELSE
         Delay_of_e_injection_T_cntr_right  = 1.0e-9 * Delay_of_e_injection_ns_right / delta_t_s                               ! Number of time steps
      END IF

   END IF

! for the right wall ions ----------------------------
   if (iBeamInjectFlag_right.eq.0) then
      if (Rank_of_process.eq.0) print '(/2x,"The ion injection at the right wall is turned off ...")'
      inject_i_every_this_many_timesteps_right = 1

   else
      temp = (iBeam_J_Am2_right / (Qs(2) * e_Cl)) / (N_plasma_m3 * delta_x_m / (delta_t_s * dble(N_of_particles_cell)))
      N_i_to_inject_total_right = int(temp)

      if (N_i_to_inject_total_right.gt.0) then
         const_N_i_to_inject_by_proc_right = N_i_to_inject_total_right / (N_of_processes-1)     ! each node will emit at least this number of ions
         variable_N_i_to_inject_right = mod(N_i_to_inject_total_right, N_of_processes-1)        ! leftover, will be distributed between all client nodes
         inject_i_every_this_many_timesteps_right = 1

      else if (temp.ge.0.001_8) then  ! minimal permitted current is one macroparticle injected every 1000 timesteps
         const_N_i_to_inject_by_proc_right = 0                               ! only one macroparticle will be injected by one (each time different) process
         variable_N_i_to_inject_right = 1                                    ! every inject_every_this_many_timesteps_right timesteps
         inject_i_every_this_many_timesteps_right = int(1.0_8 / temp) + 1    ! note that skip_timesteps_beam_inject_right >= 2

      else
         const_N_i_to_inject_by_proc_right = 0          ! if the beam was requested but its current is too small
         variable_N_i_to_inject_right = 0               ! nothing will be injected
         inject_i_every_this_many_timesteps_right = 1   !

      end if

      if (Rank_of_process.gt.0) then
         call mpi_barrier(MPI_COMM_WORLD, ierr)
         print '(2x,"Process ",i3," will inject at least ",i5," ion macroparticles at the RIGHT boundary EACH timestep")', &
         & Rank_of_process, const_N_e_to_inject_by_proc_right
      else
         ! for the server process print out the total number of macroparticles to be injected
         print '(/2x,"Constant REQUESTED flux of ion injection from the RIGHT WALL is ",e14.7," m^-2 s^-1")', (iBeam_J_Am2_right / e_Cl)
         print '(2x,"TOTAL INTEGER number of ion macroparticles to be injected from the RIGHT boundary at each timestep is : ",i7)', N_i_to_inject_total_right
         print '(2x,"TOTAL REAL    number of ion macroparticles to be injected from the RIGHT boundary at each timestep is : ",e12.5)', temp
         print '(2x,"COMMON injection leftover distributed between all nodes is : ",i5)', variable_N_i_to_inject_right
         print '(2x,"The leftover macroparticle(s) will be injected every ",i5," timestep(s)")', inject_i_every_this_many_timesteps_right
         print '("-------")'
         call mpi_barrier(MPI_COMM_WORLD, ierr)   ! place the synchronization point here to print server's output before the clients'
      end if

      VX_i_beam_right = sqrt(iBeam_energy_eV_right / (T_e_eV * Ms(2))) * N_box_vel

      ! if the delay of the start of injection is negative then it is given in timesteps, assumed to be integer !!!
      if (Delay_of_i_injection_ns_right.lt.0.0_8) then
         if (Rank_of_process.eq.0) print &
         &'(/2x,"The delay of ion injection is given in timesteps (negative value of Delay_of_i_injection_ns_right was read from the file)")'
         Delay_of_i_injection_T_cntr_right = - Delay_of_i_injection_ns_right
         Delay_of_i_injection_ns_right = Delay_of_i_injection_T_cntr_right * delta_t_s * 1.0e9

         if (Rank_of_process.eq.0) print '(/2x,"The recalculated delay of ion injection is : ", f10.3)', Delay_of_i_injection_ns_right

      else
         Delay_of_i_injection_T_cntr_right  = 1.0d-9 * Delay_of_i_injection_ns_right / delta_t_s                               ! Number of time steps
      end if

   end if

END SUBROUTINE INITIATE_INJECTION

!-------------------------------------------------------------------------------------------------
!> the subroutine injects electrons with the constant rate at both walls
!! this subroutine must be called after the final push before rearranging of the particle arays.
!! therefore it is not necessary to modify the mesh charge density here
SUBROUTINE INJECT_ELECTRONS_AT_WALLS

   USE ParallelOperationValues
   USE CurrentProblemValues, ONLY : T_cntr
   USE ElectronInjection

   IMPLICIT NONE

   INTEGER var_add, N_clients, T_cntr_flashed, k, T_inject_k

! left wall -------------------------------
   N_clients = (N_of_processes-1) * inject_e_every_this_many_timesteps_left
   T_cntr_flashed = MOD(T_cntr, N_clients) + 1  ! periodically grows from 1 to N_clients

   IF ((eBeamInjectFlag_left.NE.0).AND.(T_cntr.GE.Delay_of_e_injection_T_cntr_left)) THEN
      var_add = 0

      DO k = 0, variable_N_e_to_inject_left-1                ! this cycle works only if variable_N_e_to_inject_left>0
         ! find "timestep" suitable for injection for given client
         T_inject_k = Rank_of_process * inject_e_every_this_many_timesteps_left + k
         IF (T_inject_k.GT.N_clients) T_inject_k = T_inject_k - N_clients

         IF (T_inject_k.EQ.T_cntr_flashed) THEN
            var_add=1
            EXIT
         END IF
      END DO
      CALL INJECT_N_ELECTRONS_LEFT_WALL(const_N_e_to_inject_by_proc_left + var_add)
   END IF

! right wall -------------------------------
   N_clients = (N_of_processes-1) * inject_e_every_this_many_timesteps_right
   T_cntr_flashed = MOD(T_cntr, N_clients) + 1  ! periodically grows from 1 to N_clients

   IF ((eBeamInjectFlag_right.NE.0).AND.(T_cntr.GE.Delay_of_e_injection_T_cntr_right)) THEN
      var_add = 0

      DO k = 0, variable_N_e_to_inject_right-1                ! this cycle works only if variable_N_e_to_inject_right>0
         ! find "timestep" suitable for injection for given client
         T_inject_k = Rank_of_process * inject_e_every_this_many_timesteps_right + k
         IF (T_inject_k.GT.N_clients) T_inject_k = T_inject_k - N_clients

         IF (T_inject_k.EQ.T_cntr_flashed) THEN
            var_add=1
            EXIT
         END IF
      END DO
      CALL INJECT_N_ELECTRONS_RIGHT_WALL(const_N_e_to_inject_by_proc_right + var_add)
   END IF

END SUBROUTINE INJECT_ELECTRONS_AT_WALLS

!-------------------------------------------------------------------------------------------------
!> the subroutine injects electrons with the constant rate at both walls
!! this subroutine must be called after the final push before rearranging of the particle arays.
!! therefore it is not necessary to modify the mesh charge density here
SUBROUTINE inject_ions_at_walls

   USE ParallelOperationValues
   USE CurrentProblemValues, ONLY : T_cntr
   USE ElectronInjection

   IMPLICIT NONE

   INTEGER var_add, N_clients, T_cntr_flashed, k, T_inject_k

! left wall -------------------------------
   N_clients = (N_of_processes-1) * inject_i_every_this_many_timesteps_left
   T_cntr_flashed = MOD(T_cntr, N_clients) + 1  ! periodically grows from 1 to N_clients

   IF ((iBeamInjectFlag_left.NE.0).AND.(T_cntr.GE.Delay_of_i_injection_T_cntr_left)) THEN
      var_add = 0

      DO k = 0, variable_N_i_to_inject_left-1                ! this cycle works only if variable_N_i_to_inject_left>0
         ! find "timestep" suitable for injection for given client
         T_inject_k = Rank_of_process * inject_i_every_this_many_timesteps_left + k
         IF (T_inject_k.GT.N_clients) T_inject_k = T_inject_k - N_clients

         IF (T_inject_k.EQ.T_cntr_flashed) THEN
            var_add=1
            EXIT
         END IF
      END DO
      CALL inject_N_ions_left_wall(const_N_i_to_inject_by_proc_left + var_add)
   END IF

! right wall -------------------------------
   N_clients = (N_of_processes-1) * inject_e_every_this_many_timesteps_right
   T_cntr_flashed = MOD(T_cntr, N_clients) + 1  ! periodically grows from 1 to N_clients

   IF ((iBeamInjectFlag_right.NE.0).AND.(T_cntr.GE.Delay_of_i_injection_T_cntr_right)) THEN
      var_add = 0

      DO k = 0, variable_N_i_to_inject_right-1                ! this cycle works only if variable_N_i_to_inject_right>0
         ! find "timestep" suitable for injection for given client
         T_inject_k = Rank_of_process * inject_i_every_this_many_timesteps_right + k
         IF (T_inject_k.GT.N_clients) T_inject_k = T_inject_k - N_clients

         IF (T_inject_k.EQ.T_cntr_flashed) THEN
            var_add=1
            EXIT
         END IF
      END DO
      CALL inject_N_ions_right_wall(const_N_i_to_inject_by_proc_right + var_add)
   END IF

END SUBROUTINE inject_ions_at_walls

!-------------------------------------------------------------------------------------------------
! the subroutine injects electrons with the constant rate at the left wall
! this subroutine must be called after the final push before rearranging of the particle arays.
! therefore it is not necessary to modify the mesh charge density here
!
SUBROUTINE INJECT_N_ELECTRONS_LEFT_WALL(N_to_inject)

   USE ParallelOperationValues
   USE CurrentProblemValues
   USE ElectronInjection
   USE Diagnostics, ONLY : Rate_energy_leftemit, Rate_number_leftemit, NVX_mesh

   IMPLICIT NONE

   INTEGER N_to_inject

   REAL(8) x             ! coordinate of refluxed particle (after refluxing)
   REAL(8) vx, vy, vz    ! velocity components of refluxed particle
   REAL(8) v2            ! absolute value of velocity (squared)

   INTEGER ALLOC_ERR, j
   INTEGER left_node, right_node
   REAL(8) RAN

! left wall -------------------------------

   DO j = 1, N_to_inject

      Q_left = Q_left - Qs(1)

! get the new velocities for electrons -----------------
      vy = 0.0_8
      vz = 0.0_8
      IF (eBeamInjectFlag_left.EQ.1) THEN
         vx = VX_e_beam_left                        ! cold beam
! calculate the new coordinate
!
! note, if only one particle must be injected by all processes
! the coordinate of this particle will be jumping depending on the rank of the processes
! which is bad
!
!  basically, at this point the cold beam is not working properly #############
!
!        x = (Rank_of_process + (N_of_processes - 1) * (j - 1) - 0.5_8) * Gap_betw_e_part_left
!
         x = RAN(I_random_seed) * vx * KVx
      ELSE
         CALL GetInjMaxwellVelocity(vx)        ! warm beam, get the new x-velocity according to the maxwell distribution
         vx = vx * SQRT(eBeam_energy_eV_left / T_e_eV)    ! renormalize velocity according to the beam's temperature
! calculate the new coordinate
         x = 0.0_8
!####        x = RAN(I_random_seed) * vx * KVx
      END IF
! calculate the squared absolute velocity
      v2 = vx * vx + vy * vy + vz * vz
! Save the parameters of the refluxed electron in the linked list
      IF (ASSOCIATED(Current_electron)) THEN
         N_inject(1)            = N_inject(1) + 1
         Rate_number_leftemit(1) = Rate_number_leftemit(1) + 1             !
         NVX_mesh(0, 1) = NVX_mesh(0, 1) + 1.
         Rate_energy_leftemit(1) = Rate_energy_leftemit(1) + v2            !
         Current_electron%X       = x
         Current_electron%VX      = vx
         Current_electron%VY      = vy
         Current_electron%VZ      = vz
         Current_electron%AX      = 0.0_8
         Current_electron%Tag     = eTag_Emit_Left                         !
         CALL ADD_EMITTED_E_TO_LEFT_DF(vx, vy, vz)
         ALLOCATE(Current_electron%next, STAT = ALLOC_ERR)
         IF (ALLOC_ERR.NE.0) THEN
            PRINT '(/2x,"Process ",i3," : Error in INJECT_N_ELECTRONS_LEFT_WALL:")', Rank_of_process
            PRINT  '(2x,"Error in ALLOCATE Current_electron%next !!!")'
            PRINT  '(2x,"Program will be terminated now :(")'
            STOP
         END IF
         Current_electron => Current_electron%next
         NULLIFY(Current_electron%next)
      ELSE
         PRINT '(/2x,"Process ",i3," : Error in INJECT_N_ELECTRONS_LEFT_WALL:")', Rank_of_process
         PRINT  '(2x,"Current_electron is NOT associated!")'
         PRINT  '(2x,"The program will be terminated now :(")'
         STOP
      END IF

!       left_node  = INT(x)
!       right_node = left_node + 1
!       Q_strm_spec(left_node, 1)  = Q_strm_spec(left_node, 1)  + right_node - x    ! s = 1
!       Q_strm_spec(right_node, 1) = Q_strm_spec(right_node, 1) + x - left_node     ! s = 1
      Q_strm_spec(0, 1) = Q_strm_spec(0, 1)  + 1.0_8

   END DO

END SUBROUTINE INJECT_N_ELECTRONS_LEFT_WALL

!-------------------------------------------------------------------------------------------------
! the subroutine injects electrons with the constant rate at the right wall
! this subroutine must be called after the final push before rearranging of the particle arays.
! therefore it is not necessary to modify the mesh charge density here
!
SUBROUTINE INJECT_N_ELECTRONS_RIGHT_WALL(N_to_inject)

   USE ParallelOperationValues
   USE CurrentProblemValues
   USE ElectronInjection
   USE Diagnostics, ONLY : Rate_energy_rightemit, Rate_number_rightemit, NVX_mesh

   IMPLICIT NONE

   INTEGER N_to_inject

   REAL(8) x             ! coordinate of refluxed particle (after refluxing)
   REAL(8) vx, vy, vz    ! velocity components of refluxed particle
   REAL(8) v2            ! absolute value of velocity (squared)

   INTEGER ALLOC_ERR, j
   INTEGER left_node, right_node
   REAL(8) RAN

! right wall -------------------------------

   DO j = 1, N_to_inject

      Q_right =  Q_right - Qs(1)

! get the new velocities for electrons -----------------
      vy = 0.0_8
      vz = 0.0_8
      IF (eBeamInjectFlag_right.EQ.1) THEN
         vx = - VX_e_beam_right                        ! cold beam
! calculate the new coordinate
!
! note, if only one particle must be injected by all processes
! the coordinate of this particle will be jumping depending on the rank of the processes
! which is bad
!
!  basically, at this point the cold beam is not working properly #############
!
!        x = N_cells - (Rank_of_process + (N_of_processes - 1) * (j - 1) - 0.5_8) * Gap_betw_part_right
         x = N_cells + RAN(I_random_seed) * vx * KVx
      ELSE
         CALL GetInjMaxwellVelocity(vx)        ! warm beam, get the new x-velocity according to the maxwell distribution
         vx = - vx * SQRT(eBeam_energy_eV_right / T_e_eV)    ! renormalize velocity according to the beam's temperature
! calculate the new coordinate
         x = DBLE(N_cells)
!####        x = N_cells + RAN(I_random_seed) * vx * KVx
      END IF
! calculate the squared absolute velocity
      v2 = vx * vx + vy * vy + vz * vz
! Save the parameters of the refluxed electron in the linked list
      IF (ASSOCIATED(Current_electron)) THEN
         N_inject(1)            = N_inject(1) + 1
         Rate_number_rightemit(1) = Rate_number_rightemit(1) + 1             !
         NVX_mesh(N_cells, 1) = NVX_mesh(N_cells, 1) - 1.
         Rate_energy_rightemit(1) = Rate_energy_rightemit(1) + v2            !
         Current_electron%X       = x
         Current_electron%VX      = vx
         Current_electron%VY      = vy
         Current_electron%VZ      = vz
         Current_electron%AX      = 0.0_8
         Current_electron%Tag     = eTag_Emit_Right                         !
         CALL ADD_EMITTED_E_TO_RIGHT_DF(vx, vy, vz)
         ALLOCATE(Current_electron%next, STAT = ALLOC_ERR)
         IF (ALLOC_ERR.NE.0) THEN
            PRINT '(/2x,"Process ",i3," : Error in INJECT_N_ELECTRONS_RIGHT_WALL:")', Rank_of_process
            PRINT  '(2x,"Error in ALLOCATE Current_electron%next !!!")'
            PRINT  '(2x,"Program will be terminated now :(")'
            STOP
         END IF
         Current_electron => Current_electron%next
         NULLIFY(Current_electron%next)
      ELSE
         PRINT '(/2x,"Process ",i3," : Error in INJECT_N_ELECTRONS_RIGHT_WALL:")', Rank_of_process
         PRINT  '(2x,"Current_electron is NOT associated!")'
         PRINT  '(2x,"The program will be terminated now :(")'
         STOP
      END IF

!! account for the density of emitted particle
!     left_node  = MIN(INT(x), N_cells-1)
!     right_node = left_node + 1
!     Q_strm_spec(left_node, 1)  = Q_strm_spec(left_node, 1)  + right_node - x    ! s = 1
!     Q_strm_spec(right_node, 1) = Q_strm_spec(right_node, 1) + x - left_node     ! s = 1

      Q_strm_spec(N_cells, 1) = Q_strm_spec(N_cells, 1) + 1.0_8

   END DO

END SUBROUTINE INJECT_N_ELECTRONS_RIGHT_WALL

!-------------------------------------------------------------------------------------------------
! Rewrite the inject_n_electrons_left_wall and inject_n_electrons_right_wall subroutines for
! ions, and make sure to use lower case for all the commands
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
! the subroutine injects ions with the constant rate at the left wall
! this subroutine must be called after the final push before rearranging of the particle arays.
! therefore it is not necessary to modify the mesh charge density here
subroutine inject_N_ions_left_wall(N_to_inject)

   use ParallelOperationValues
   use CurrentProblemValues
   use ElectronInjection
   use Diagnostics, only : Rate_energy_leftemit, Rate_number_leftemit, NVX_mesh

   implicit none

   integer N_to_inject

   real(8) x             ! coordinate of refluxed particle (after refluxing)
   real(8) vx, vy, vz    ! velocity components of refluxed particle
   real(8) v2            ! absolute value of velocity (squared)

   integer ALLOC_ERR, j
   integer left_node, right_node
   real(8) RAN

! left wall -------------------------------

   do j = 1, N_to_inject

      Q_left = Q_left - Qs(2)

! get the new velocities for electrons -----------------
      vy = 0.0_8
      vz = 0.0_8
      if (iBeamInjectFlag_left.eq.1) then
         vx = VX_i_beam_left                        ! cold beam
! calculate the new coordinate
!
! note, if only one particle must be injected by all processes
! the coordinate of this particle will be jumping depending on the rank of the processes
! which is bad
!
!  basically, at this point the cold beam is not working properly #############
!
!        x = (Rank_of_process + (N_of_processes - 1) * (j - 1) - 0.5_8) * Gap_betw_e_part_left
!
         x = RAN(I_random_seed) * vx * KVx
      else
         call GetInjMaxwellVelocity(vx)        ! warm beam, get the new x-velocity according to the maxwell distribution
         vx = vx * sqrt(iBeam_energy_eV_left / (Ms(2) * T_e_eV))    ! renormalize velocity according to the beam's temperature
! calculate the new coordinate
         x = 0.0_8
!####        x = RAN(I_random_seed) * vx * KVx
      end if
! calculate the squared absolute velocity
      v2 = vx * vx + vy * vy + vz * vz
! Save the parameters of the refluxed electron in the linked list
      if (associated(Current_ion)) then
         N_inject(2)            = N_inject(2) + 1
         Rate_number_leftemit(2) = Rate_number_leftemit(2) + 1             !
         NVX_mesh(0, 2) = NVX_mesh(0, 2) + 1.
         Rate_energy_leftemit(2) = Rate_energy_leftemit(2) + v2            !
         Current_ion%X       = x
         Current_ion%VX      = vx
         Current_ion%VY      = vy
         Current_ion%VZ      = vz
         Current_ion%AX      = 0.0_8
         Current_ion%Tag     = 0                         !
         allocate(Current_ion%next, stat = ALLOC_ERR)
         if (ALLOC_ERR.NE.0) then
            print '(/2x,"Process ",i3," : Error in inject_N_ions_left_wall:")', Rank_of_process
            print  '(2x,"Error in ALLOCATE Current_ion%next !!!")'
            print  '(2x,"Program will be terminated now :(")'
            stop
         end if
         Current_ion => Current_ion%next
         nullify(Current_ion%next)
      else
         print '(/2x,"Process ",i3," : Error in inject_N_ions_left_wall:")', Rank_of_process
         print  '(2x,"Current_ion is NOT associated!")'
         print  '(2x,"The program will be terminated now :(")'
         stop
      end if

      Q_strm_spec(0, 2) = Q_strm_spec(0, 2)  + 1.0_8

   end do

end subroutine inject_N_ions_left_wall

!-------------------------------------------------------------------------------------------------
! the subroutine injects ions with the constant rate at the right wall
! this subroutine must be called after the final push before rearranging of the particle arays.
! therefore it is not necessary to modify the mesh charge density here
subroutine inject_N_ions_right_wall(N_to_inject)

   use ParallelOperationValues
   use CurrentProblemValues
   use ElectronInjection
   use Diagnostics, only : Rate_energy_rightemit, Rate_number_rightemit, NVX_mesh

   implicit none

   integer N_to_inject

   real(8) x             ! coordinate of refluxed particle (after refluxing)
   real(8) vx, vy, vz    ! velocity components of refluxed particle
   real(8) v2            ! absolute value of velocity (squared)

   integer ALLOC_ERR, j
   integer left_node, right_node
   real(8) RAN

! right wall -------------------------------

   do j = 1, N_to_inject

      Q_right =  Q_right - Qs(2)

! get the new velocities for electrons -----------------
      vy = 0.0_8
      vz = 0.0_8
      if (iBeamInjectFlag_right.eq.1) then
         vx = - VX_i_beam_right                        ! cold beam
! calculate the new coordinate
!
! note, if only one particle must be injected by all processes
! the coordinate of this particle will be jumping depending on the rank of the processes
! which is bad
!
!  basically, at this point the cold beam is not working properly #############
!
!        x = N_cells - (Rank_of_process + (N_of_processes - 1) * (j - 1) - 0.5_8) * Gap_betw_part_right
         x = N_cells + RAN(I_random_seed) * vx * KVx
      else
         call GetInjMaxwellVelocity(vx)        ! warm beam, get the new x-velocity according to the maxwell distribution
         vx = - vx * sqrt(iBeam_energy_eV_right / (Ms(2) * T_e_eV))    ! renormalize velocity according to the beam's temperature
         ! calculate the new coordinate
         x = dble(N_cells)
!####        x = N_cells + RAN(I_random_seed) * vx * KVx
      end if
      ! calculate the squared absolute velocity
      v2 = vx * vx + vy * vy + vz * vz
      ! Save the parameters of the emitted ion in the linked list
      if (associated(Current_ion)) then
         N_inject(2)            = N_inject(2) + 1
         Rate_number_rightemit(2) = Rate_number_rightemit(2) + 1             !
         NVX_mesh(N_cells, 2) = NVX_mesh(N_cells, 2) - 1.
         Rate_energy_rightemit(2) = Rate_energy_rightemit(2) + v2            !
         Current_ion%X       = x
         Current_ion%VX      = vx
         Current_ion%VY      = vy
         Current_ion%VZ      = vz
         Current_ion%AX      = 0.0_8
         Current_ion%Tag     = 0                         !
         allocate(Current_ion%next, stat = ALLOC_ERR)
         if (ALLOC_ERR.NE.0) then
            print '(/2x,"Process ",i3," : Error in inject_N_ions_right_wall:")', Rank_of_process
            print  '(2x,"Error in ALLOCATE Current_ion%next !!!")'
            print  '(2x,"Program will be terminated now :(")'
            STOP
         end if
         Current_ion => Current_ion%next
         nullify(Current_ion%next)
      else
         print '(/2x,"Process ",i3," : Error in inject_N_ions_right_wall:")', Rank_of_process
         print  '(2x,"Current_ion is NOT associated!")'
         print  '(2x,"The program will be terminated now :(")'
         stop
      end if

      Q_strm_spec(N_cells, 1) = Q_strm_spec(N_cells, 1) + 1.0_8

   end do

end subroutine inject_N_ions_right_wall