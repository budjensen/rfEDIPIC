!----------------------------------
module RF_diagnostics

   integer, allocatable :: RF_Diagnostic_steps(:)  !! Timesteps RF diagnostics are printed
   integer current_RF_step                         !! Number of the current RF diagnostic step (starts at 0)

end module RF_diagnostics

!----------------------------------
subroutine initialize_RF_diagnostics()
   !! This subroutine initializes the RF diagnostics by setting the current RF step to 1
   use CurrentProblemValues
   use RF_diagnostics
   use heating_diagnostics
   implicit none
   CHARACTER (77) buf
   INTEGER :: exists, i, ierr

! read / write the data file
   inquire (file = 'ssc_rfdiagnostics.dat', EXIST = exists)
   call MPI_BARRIER(MPI_COMM_WORLD, ierr)

   open (9, FILE = 'ssc_diagnostics.dat')
   if (exists) then

      if (Rank_of_process.EQ.0) then
         PRINT '(2x,"Process ",i3," : RF Diagnostics control data file is found. Reading the data file...")', Rank_of_process
      END IF

      read (9, '(A77)') buf ! '("======================== RF Period Snapshot Control =========================")')
      read (9, '(A77)') buf ! '("-------d----- RF Diagnostics Control Switch ( 1/0 = on/off ) ----------------")')
      read (9, '(7x,i1)') Flag_RF_diagnostics
      read (9, '(A77)') buf ! '("======================= RF Snapshot Averaging Control =======================")')
      read (9, '(A77)') buf ! '("--dddddd----- RF Period of first collection ---------------------------------")')
      read (9, '(2x,i6)') Start_RF_collection
      read (9, '(A77)') buf ! '("---ddddd----- Number of RF periods to average -------------------------------")')
      read (9, '(3x,i5)') N_RF_periods_to_avg
      read (9, '(A77)') buf ! '("------dd----- Number of collection times within one RF period ---------------")')
      read (9, '(6x,i2)') N_collection_in_RF_period
      read (9, '(A77)') buf ! '("---0.dddddd-- Collection times, in radians ( 0 <= t < 2pi = 6.283185 ) ------")')
      do i = 1, N_collection_in_RF_period
         read (9, '(3x,f8.6)') RF_collection_time
         if ((RF_collection_time.ge.0.0_8).and.(RF_collection_time.lt.6.283185_8)) then
            RF_collection_times(N_of_probes) = RF_collection_time
         else
            print '(/2x,"Error in initialize_RF_diagnostics: Collection time must be in range [0, 2pi=6.283185] !!!")'
            print  '(2x,"Program will be terminated now :(")'
            STOP
         end if
      END DO

   else

!======================= Time dependencies control =======================
      WriteOut_step            = 100
      flag_Out                 = 1
      WriteStart_step          = 0
      flag_Start             = 1
      WriteAvg_step            = 1
      flag_Avg                 = 1
      TextOut_avg_prds_to_skip = 0
      N_of_probes_given        = 0
!======================= Snapshot control =======================
      N_of_snap_groups    = 1
      Aprx_snap_start_ns  = 0.0_8
      Aprx_snap_finish_ns = 1.0_8
      Aprx_n_of_snaps     = 2
      Flag_pp_creation    = 0
      N_to_skip           = 100
      N_to_skip_left      = 0
      N_to_skip_right     = 0
      N_to_skip_ion       = 100
      Ve_x_max            = 3.0_8
      Ve_yz_max           = 3.0_8
      N_of_breakpoints    = 0
      Ee_max_eV           = 100.0_8
      Ei_max_eV           = 20.0_8
      N_of_E_breakpoints  = 0
      Ei_wall_max_eV      = 100.0_8
      N_E_bins            = 50
      flag_df_temp(1,1)   = 1
      flag_df_temp(1,2)   = 1
      flag_df_temp(1,3)   = 1
      flag_df_temp(1,4)   = 1
      flag_df_temp(2,1)   = 1
      flag_df_temp(2,2)   = 1
      flag_df_temp(2,3)   = 1
      flag_df_temp(2,4)   = 1

      IF (Rank_of_process.EQ.0) THEN

         PRINT '(/2x,"Process ",i3," : File ssc_diagnostics.dat not found. Use the default settings ...")', Rank_of_process
         PRINT '(2x,"Process ",i3," : Create ssc_diagnostics.dat file . . .")', Rank_of_process

         WRITE (9, '("********************* TIME DEPENDENCIES CREATION CONTROL ********************")')
         WRITE (9, '("------------- Diagnostic output interval, t_flag ----------------------------")')
         WRITE (9, '("--            If t_flag <= 0:   t_output * 10^(t_flag) in RF periods --------")')
         WRITE (9, '("--            If t_flag  > 1:   t_output > 0 given in timesteps -------------")')
         WRITE (9, '("-dddddddd-#d-                   t_output < 0 given in plasma periods --------")')
         WRITE (9, '(1x,i7,1x,i2)') WriteOut_step, flag_Out
         WRITE (9, '("------------- Diagnostic collection start, t_flag ---------------------------")')
         WRITE (9, '("--            If t_flag <= 0:   t_start * 10^(t_flag) in RF periods ---------")')
         WRITE (9, '("--            If t_flag  > 1:   t_start > 0 given in timesteps --------------")')
         WRITE (9, '("-dddddddd-#d-                   t_start < 0 given in plasma periods ---------")')
         WRITE (9, '(1x,i8,1x,i2)') WriteStart_step, flag_Start
         WRITE (9, '("------------- Diagnostic average window, t_flag - (NOTE: t_avg <= t_output) -")')
         WRITE (9, '("--            If t_flag <= 0:   t_avg * 10^(t_flag) in RF periods -----------")')
         WRITE (9, '("--            If t_flag  > 1:   t_avg > 0 given in timesteps ----------------")')
         WRITE (9, '("-dddddddd-#d-                   t_avg < 0 given in plasma periods -----------")')
         WRITE (9, '(1x,i7,1x,i2)') WriteAvg_step, flag_Avg
         WRITE (9, '("--dddddd----- Skip periods of averaging between text outputs (>=0) ----------")')
         WRITE (9, '(2x,i6)') TextOut_avg_prds_to_skip
         WRITE (9, '("-----ddd----- Number of probes ( no probes if <= 0 ) ------------------------")')
         WRITE (9, '(5x,i3)') N_of_probes_given
         WRITE (9, '("--dddddd.ddd- Probe coordinates (node number if>0 or millimeters if<0) ------")')
         WRITE (9, '("************************ SNAPSHOTS CREATION CONTROL *************************")')
         WRITE (9, '("------dd----- Number of groups of snapshots ( >= 0 )-------------------------")')
         WRITE (9, '(6x,i2)') N_of_snap_groups
         WRITE (9, '("----------- start (ns) ---------- finish (ns) ---------- number -------------")')
         WRITE (9, '("------------dddddd.ddd------------dddddd.ddd--------------dddd---------------")')
         WRITE (9, '(12x,f10.3,12x,f10.3,14x,i4)') Aprx_snap_start_ns, Aprx_snap_finish_ns, Aprx_n_of_snaps
         WRITE (9, '("********************** PHASE PLANES CREATION CONTROL ************************")')
         WRITE (9, '("-------d----- Create phase planes ? (1/0) -----------------------------------")')
         WRITE (9, '(7x,i1)') Flag_pp_creation
         WRITE (9, '("-ddd-ddd-ddd-ddd--- Number of particles to skip (bulk/left/right/ion, >=0) --")')
         WRITE (9, '(4(1x,i3))') N_to_skip, N_to_skip_left, N_to_skip_right, N_to_skip_ion
         WRITE (9, '("************* VELOCITY DISTRIBUTION FUNCTIONS CREATION CONTROL **************")')
         WRITE (9, '("--dddddd.ddd- Maximal x-velocity for e-v_x-distribution (in V_therm_e) ------")')
         WRITE (9, '(2x,f10.3)') Ve_x_max
         WRITE (9, '("--dddddd.ddd- Maximal y,z-velocity for e-v_y,z-distribution (in V_therm_e) --")')
         WRITE (9, '(2x,f10.3)') Ve_yz_max
         WRITE (9, '("------dd----- Number of breakpoints (no evdfs if<0, all system used if 0) ---")')
         WRITE (9, '(6x,i2)') N_of_breakpoints
         WRITE (9, '("--dddddd.ddd- Breakpoints (ascend., node number if>0 or millimeters if <0) --")')
         WRITE (9, '("************** ENERGY DISTRIBUTION FUNCTIONS CREATION CONTROL ***************")')
         WRITE (9, '("--dddddd.ddd- Maximal energy for electron e-distribution (in eV) ------------")')
         WRITE (9, '(2x,f10.3)') Ee_max_eV
         WRITE (9, '("--dddddd.ddd- Maximal energy for ion e-distribution (in eV) -----------------")')
         WRITE (9, '(2x,f10.3)') Ei_max_eV
         WRITE (9, '("------dd----- Number of breakpoints (no edfs if<0, all system used if 0) ----")')
         WRITE (9, '(6x,i2)') N_of_E_breakpoints
         WRITE (9, '("--dddddd.ddd- Breakpoints (ascend., node number if>0 or millimeters if <0) --")')
         WRITE (9, '("--dddddd.ddd- Maximal energy for ion wall e-distribution (in eV) ------------")')
         WRITE (9, '(2x,f10.3)') Ei_wall_max_eV
         WRITE (9, '("-----ddd----- Number of bins for e-distributions (>0) -----------------------")')
         WRITE (9, '(5x,i3)') N_E_bins
         WRITE (9, '("********************** DISTRIBUTION FUNCTION SWITCHES ***********************")')
         WRITE (9, '("--d--d--d--d- Velocity ( e-vx, e-vy, e-vz, i-vx, | 1/0 = on/off ) -----------")')
         WRITE (9, '(4(2x,i1))') flag_df_temp(1,1), flag_df_temp(1,2), flag_df_temp(1,3), flag_df_temp(1,4)
         WRITE (9, '("--d--d--d--d- Energy ( e, i, i l-wall, i r-wall | 1/0 = on/off ) ------------")')
         WRITE (9, '(4(2x,i1))') flag_df_temp(1,1), flag_df_temp(1,2), flag_df_temp(1,3), flag_df_temp(1,4)

      END IF

   END IF

   CLOSE (9, STATUS = 'KEEP')

   ! Set the current RF step to 0
   current_RF_step = 0

   ! Calculate the RF diagnostic steps
   call calc_RF_diagnostic_steps()


end subroutine initialize_RF_diagnostics

subroutine calc_RF_diagnostic_steps()
   !! This subroutine calculates the time steps at which the RF diagnostics will be written
   use CurrentProblemValues, only : delta_t_s, f_rf_Hz, t_start_s
   use RF_diagnostics, only: RF_Diagnostic_steps
   implicit none
   integer RF_period_start

   ! Calculate the timestep where the rf source turns on
   RF_period_start = int(t_start_s / delta_t_s) + 1



end subroutine calc_RF_diagnostic_steps

subroutine allocate_RF_diag_arrays()
   use RF_diagnostics, only : RF_Diagnostic_steps
   implicit none
   integer :: ALLOC_ERR

   ! Allocate the heating array
   if (.not.allocated(RF_Diagnostic_steps)) then
      allocate(RF_Diagnostic_steps(), stat=ALLOC_ERR)
      if(ALLOC_ERR.ne.0) then
         print *, 'Error in ALLOCATE RF_Diagnostic_steps !!!'
         print *, 'The program will be terminated now :('
         stop
      end if
   end if

end subroutine allocate_RF_diag_arrays

subroutine write_RF_diagnostics()
   !! This subroutine prints out the instantaneous system potential, electric field,
   !! displacement current, and electron and ion densities, fluxes, and currents at
   !! each node at the time steps specified by RF_Diagnostic_steps.
   use mpi
   use ParallelOperationValues
   use CurrentProblemValues
   use Diagnostics, only : flux_m2s1
   use RF_diagnostics
   use heating_diagnostics
   implicit none

   character(23) rf_diag_filename
   character(4) rf_diag_number_txt
   integer filename_blnks

   real(8) Ne_m3, Ni_m3             ! These are multiplied by 2 at the first and last node since N_scl_m3 = N_plasma_m3 / N_of_particles_cell.
                                    ! To tally up the total number of particles correctly we need to ensure each node point counts
                                    ! the particles in one cell, except for the end nodes which count the particles in half a cell
                                    ! (thus the number of particles is correctly tallied)
   real(8) Jx_Am2, Jy_Am2, Jz_Am2   ! These are multiplied by 2 at the first and last node since N_scl_m3 = N_plasma_m3 / N_of_particles_cell.
                                    ! To tally up the total number of particles correctly we need to ensure each node point counts
                                    ! the particles in one cell, except for the end nodes which count the particles in half a cell
                                    ! (thus the number of particles is correctly tallied)

   integer ierr, ALLOC_ERR, DEALLOC_ERR
   integer, allocatable :: ibufer(:)      !! send/receive buffer for integer data
   real(8), allocatable :: dbufer(:)      !! send/receive buffer for real(8) data
   integer, allocatable :: ibufer2(:)     !! send/receive buffer for integer data
   real(8), allocatable :: dbufer2(:)     !! send/receive buffer for real(8) data

   if (T_cntr.ne.RF_Diagnostic_steps(current_RF_step + 1)) return

   ! Increment the current RF step
   current_RF_step = current_RF_step + 1

   ! Allocate buffers for send/receive
   IF (.NOT.ALLOCATED(ibufer)) THEN
      ALLOCATE(ibufer(1:N_nodes), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : write_RF_diagnostics : Error in ALLOCATE ibufer !!!")', Rank_of_process
         PRINT '(2x,"The program will be terminated now :(")'
         STOP
      END IF
   END IF
   IF (.NOT.ALLOCATED(ibufer2)) THEN
      ALLOCATE(ibufer2(1:N_nodes), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : write_RF_diagnostics : Error in ALLOCATE ibufer2 !!!")', Rank_of_process
         PRINT '(2x,"The program will be terminated now :(")'
         STOP
      END IF
   END IF
   IF (.NOT.ALLOCATED(rbufer)) THEN
      ALLOCATE(rbufer(1:N_nodes), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : write_RF_diagnostics : Error in ALLOCATE rbufer !!!")', Rank_of_process
         PRINT '(2x,"The program will be terminated now :(")'
         STOP
      END IF
   END IF
   IF (.NOT.ALLOCATED(rbufer2)) THEN
      ALLOCATE(rbufer2(1:N_nodes), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : write_RF_diagnostics : Error in ALLOCATE rbufer2 !!!")', Rank_of_process
         PRINT '(2x,"The program will be terminated now :(")'
         STOP
      END IF
   END IF

   ! Send from client to server
   if (Rank_of_process.NE.0) then                                                                 ! the client >>>
      ! SPATIAL PROFILES
      ! --------------------------------
      ! Q_strm_spec(0:N_cells,1:N_spec) is already on the server
      ! F(0:N_cells) is already on the server
      ! EX(0:N_cells) is already on the server
      ! joule_heat_Wm3
      DO s = 1, N_spec
         dbufer(1:N_cells)  = joule_heat_Wm3(1:N_cells, s)
         dbufer2(1:N_cells) = 0.0_8
         CALL MPI_REDUCE(dbufer, dbufer2, N_cells, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO
      ! flux_m2s1
      do s =  1, N_spec
         dbufer(1:N_cells)  = flux_m2s1(1:N_cells, s)
         dbufer2(1:N_cells) = 0.0_8
         CALL MPI_REDUCE(dbufer, dbufer2, N_cells, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      end do
      ! --------------------------------
   
      ! ! GLOBAL DATA (none at the moment)
      ! ! --------------------------------
      ! ! prepare integer buffer for transmission
      ! ibufer( 1) = N_part(1)
      ! ibufer( 2) = N_part(2)

      ! ! prepare the real(8) buffer for transmission
      ! dbufer( 1) = VY_recent(1)
      ! dbufer( 2) = VY_recent(2)

      ! ! transmit
      ! CALL MPI_REDUCE(ibufer, ibufer2, 2, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      ! CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

      ! CALL MPI_REDUCE(dbufer, dbufer2, 2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      ! CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

      ! quit
      RETURN

   ! Receive from client
   ELSE                                                                                           ! the server >>>
      PRINT '(/2x,"Process ",i3," : ^^^^^^^^^^^^ RF Diagnostics number ",i4," will be saved now ... ^^^^^^^^^^^^")', &
      & Rank_of_process, current_RF_step
   
      ! SPATIAL PROFILES
      ! --------------------------------
      ! receive data from clients
      ! joule_heat_Wm3
      DO s = 1, N_spec
         dbufer  = 0.0_8
         dbufer2 = 0.0_8
         CALL MPI_REDUCE(dbufer2, dbufer, N_cells, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         joule_heat_Wm3(1:N_cells, s) = dbufer(1:N_cells)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO
      ! flux_m2s1
      do s =  1, N_spec
         dbufer(1:N_cells)  = 0.0_8
         dbufer2(1:N_cells) = 0.0_8
         CALL MPI_REDUCE(dbufer2, dbufer, N_cells, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         flux_m2s1(1:N_cells, s) = dbufer(1:N_cells)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      end do
      ! --------------------------------

      ! ! GLOBAL DATA (none at the moment)
      ! ! --------------------------------
      ! ! clear the buffers
      ! ibufer = 0
      ! dbufer = 0.0_8

      ! ! receive
      ! CALL MPI_REDUCE(ibufer2, ibufer, 2, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      ! CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

      ! CALL MPI_REDUCE(dbufer2, dbufer, 2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      ! CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

      ! ! restore integer data from the buffer
      ! N_part(1)                   = ibufer( 1)
      ! N_part(2)                   = ibufer( 2)

      ! ! restore real(8) data from the buffer
      ! VY_recent(1)                 = dbufer( 1)
      ! VY_recent(2)                 = dbufer( 2)

   END IF

   !----------------------------------------------------
   ! Write the RF diagnostics to file

   ! Produce new filename numbers
   WRITE (rf_diag_number_txt, '(i4)') current_RF_step       ! Conversation of integer i to character*4 number_txt
   rf_diag_number_txt = adjustl(trip(rf_diag_number_txt))    ! Align to the left left in order to
   filename_blnks = 4 - len_trim(rf_diag_number_txt)         ! calculate the number of blanks in string number_txt;
   rf_diag_number_txt = adjustr(rf_diag_number_txt)          ! then align to the right and fulfill
   rf_diag_number_txt(1:filename_blnks) = '0000'             ! substitution of leading blanks by '0'

   ! Produce period of the new filename
   current_time_s = T_cntr * delta_t_s


   ! create the filenames
   rf_diag_filename =     '_NNNN_period_TTTTTT.dat'
   rf_diag_filename(2:5) = rf_diag_number_txt
   rf_diag_filename(14:19) = rf_period_txt

   OPEN (99, FILE = rf_diag_filename)

   ! save column description
   WRITE (99, '("# col  1 is X-coordinate [m]")')

end subroutine write_RF_diagnostics

