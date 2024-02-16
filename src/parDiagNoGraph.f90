!------------------------------
!
SUBROUTINE INITIATE_DIAGNOSTICS

   use mpi
   USE ParallelOperationValues
   USE CurrentProblemValues
   USE Diagnostics
   USE Snapshots

   IMPLICIT NONE

!  INCLUDE 'mpif.h'

   LOGICAL exists
   INTEGER ALLOC_ERR !, DEALLOC_ERR

   INTEGER i      ! set of snapshots

   INTEGER N_of_snap_groups          ! number of sets of snapshots, read from file
   REAL(8) Aprx_snap_start_ns        ! approximate start of current set of snapshots [ns], read from file
   REAL(8) Aprx_snap_finish_ns       ! approximate finish of current set of snapshots [ns], read from file
   INTEGER Aprx_n_of_snaps           ! approximate number of snapshots in current set, read from file

   INTEGER T2_old, N2_old !, T_current

   INTEGER N_of_breakpoints          ! approximate number of breakpoints - right boundaries of
   ! locations for calculation of velocity distribution functions (VDF)
   ! N_of_all_vdf_locs = N_of_breakpoints + 1
   ! if N_of_breakpoints < 0 no local distributions will be created
   ! if N_of_breakpoints = 0 distribution will be created using all particles
   INTEGER N_of_E_breakpoints        ! approximate number of breakpoints - right boundaries of
   ! locations for calculation of energy distribution functions (EDF)
   ! N_of_all_edf_locs = N_of_E_breakpoints + 1
   ! if N_of_E_breakpoints < 0 no local distributions will be created
   ! if N_of_E_breakpoints = 0 distribution will be created using all particles
   REAL(8) X_breakpoint              ! breakpoint coordinate

   INTEGER N_of_probes_given         ! number of lines in data file where probe locations are described

   INTEGER timestep(1:9999)          ! array for temporary storage of moments (timesteps) of snapshots
   INTEGER breakpoint(1:200)         ! array for temporary storage of locations (node numbers) of breakpoints ## (double the size to filter
   INTEGER breakpoint_E(1:200)       ! array for temporary storage of locations (node numbers) of breakpoints ##  input with multiple errors)
   INTEGER probe_location(1:2000)    ! array for temporary storage of locations (node numbers) of probes      ##  input with multiple errors)

   integer flag_Out              !! Flag for input file timing variable WriteOut_step
                                 !! If flag <= 0:   t_output * 10^(flag) in RF periods
                                 !! If flag  > 0:   t_output > 0 given in timesteps
                                 !!                 t_output < 0 given in plasma periods
   integer flag_Start          !! Flag for input file timing variable WriteStart_step
                                 !! If flag <= 0:   t_output * 10^(flag) in RF periods
                                 !! If flag  > 0:   t_output > 0 given in timesteps
                                 !!                 t_output < 0 given in plasma periods
   integer flag_Avg              !! Flag for input file timing variable WriteAvg_step
                                 !! If flag <= 0:   t_output * 10^(flag) in RF periods
                                 !! If flag  > 0:   t_output > 0 given in timesteps
                                 !!                 t_output < 0 given in plasma periods
   integer flag_df_temp(1:2,1:4) !! Temporary variable for holding flags for distribution function creation
                                 !! (1 [velocity], 1:4): 1 - e-vx, 2 - e-vy, 3 - e-vz, 4 - i-vx
                                 !! (2 [energy], 1:4): 1 - e, 2 - i, 3 - i l-wall, 4 - i r-wall

   real(8) steps_per_period      !! Number of timesteps per RF period
   real(8) NumStart_RFper       !! Number of RF periods in WriteStart_step
   real(8) NumOut_RFper         !! Number of RF periods in WriteOut_step
   real(8) NumAvg_RFper         !! Number of RF periods in WriteAvg_step

   CHARACTER (77) buf
   INTEGER ierr

   N_of_all_snaps = 0
   N_of_all_vdf_locs  = 0
   Flag_pp_creation = 0
   N_of_probes = 0
   T2_old = 0
   N2_old = 0

! read / write the data file
   INQUIRE (FILE = 'ssc_diagnostics.dat', EXIST = exists)
   CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

   OPEN (9, FILE = 'ssc_diagnostics.dat')
   IF (exists) THEN

      IF (Rank_of_process.EQ.0) THEN
         PRINT '(2x,"Process ",i3," : Diagnostics control data file is found. Reading the data file...")', Rank_of_process
      END IF

      READ (9, '(A77)') buf ! '("********************* TIME DEPENDENCIES CREATION CONTROL ********************")')
      READ (9, '(A77)') buf ! '("------------- Diagnostic output interval, t_flag ----------------------------")')
      READ (9, '(A77)') buf ! '("--            If t_flag <= 0:   t_output * 10^(t_flag) in RF periods --------")')
      READ (9, '(A77)') buf ! '("--            If t_flag  > 1:   t_output > 0 given in timesteps -------------")')
      READ (9, '(A77)') buf ! '("-dddddddd-#d-                   t_output < 0 given in plasma periods --------")')
      READ (9, '(1x,i8,1x,i2)') WriteOut_step, flag_Out
      READ (9, '(A77)') buf ! '("------------- Diagnostic collection start, t_flag ---------------------------")')
      READ (9, '(A77)') buf ! '("--            If t_flag <= 0:   t_start * 10^(t_flag) in RF periods ---------")')
      READ (9, '(A77)') buf ! '("--            If t_flag  > 1:   t_start > 0 given in timesteps --------------")')
      READ (9, '(A77)') buf ! '("-dddddddd-#d-                   t_start < 0 given in plasma periods ---------")')
      READ (9, '(1x,i8,1x,i2)') WriteStart_step, flag_Start
      READ (9, '(A77)') buf ! '("------------- Diagnostic average window, t_flag - (NOTE: t_avg <= t_output) -")')
      READ (9, '(A77)') buf ! '("--            If t_flag <= 0:   t_avg * 10^(t_flag) in RF periods -----------")')
      READ (9, '(A77)') buf ! '("--            If t_flag  > 1:   t_avg > 0 given in timesteps ----------------")')
      READ (9, '(A77)') buf ! '("-dddddddd-#d-                   t_avg < 0 given in plasma periods -----------")')
      READ (9, '(1x,i8,1x,i2)') WriteAvg_step, flag_Avg
      READ (9, '(A77)') buf ! '("--dddddd----- Skip periods of averaging between text outputs (>=0) ----------")')
      READ (9, '(2x,i6)') TextOut_avg_prds_to_skip
      READ (9, '(A77)') buf ! '("-----ddd----- Number of probes ( no probes if <= 0 ) ------------------------")')
      READ (9, '(5x,i3)') N_of_probes_given
      READ (9, '(A77)') buf ! '("--dddddd.ddd- Probe coordinates (node number if>0 or millimeters if<0) ------")')

      DO i = 1, N_of_probes_given
         READ (9, '(2x,f10.3)') X_breakpoint
         IF (X_breakpoint.LT.0.0_8) X_breakpoint = 0.001_8 * ABS(X_breakpoint) / delta_x_m  ! convert from millimeters
         IF ( (INT(X_breakpoint).GE.0) .AND. (INT(X_breakpoint).LE.N_cells) ) THEN ! probes must be within the domain
            N_of_probes = N_of_probes + 1
            probe_location(N_of_probes) = INT(X_breakpoint)
         END IF
      END DO   ! note that after this cycle, N_of_probes is not necessary equal to N_of_probes_given

      READ (9, '(A77)') buf ! '("************************ SNAPSHOTS CREATION CONTROL *************************")')
      READ (9, '(A77)') buf ! '("------dd----- Number of groups of snapshots ( >= 0 )-------------------------")')
      READ (9, '(6x,i2)') N_of_snap_groups
      READ (9, '(A77)') buf ! '("----------- start (ns) ---------- finish (ns) ---------- number -------------")')
      READ (9, '(A77)') buf ! '("------------dddddd.ddd------------dddddd.ddd--------------dddd---------------")')

      DO i = 1, N_of_snap_groups
         ! Read the parameters of current set of snapshot from the data file
         READ (9, '(12x,f10.3,12x,f10.3,14x,i4)') Aprx_snap_start_ns, Aprx_snap_finish_ns, Aprx_n_of_snaps    !!!
         CALL CALCULATE_SNAPSHOT_TIMES(Aprx_snap_start_ns, Aprx_snap_finish_ns, Aprx_n_of_snaps, T2_old, N2_old, timestep)
      END DO

      IF (N_of_all_snaps.GT.0) THEN ! skip further reading if no snapshots requested

         READ (9, '(A77)') buf ! ********************** PHASE PLANES CREATION CONTROL ************************")')
         READ (9, '(A77)') buf ! -------d----- Create phase planes ? (1/0) -----------------------------------")') creation of phase planes can be cancelled
         READ (9, '(7x,i1)') Flag_pp_creation
         READ (9, '(A77)') buf ! -ddd-ddd-ddd-ddd--- Number of particles to skip (bulk/left/right/ion, >=0) --")')
         READ (9, '(4(1x,i3))') N_to_skip, N_to_skip_left, N_to_skip_right, N_to_skip_ion

         READ (9, '(A77)') buf ! ************* VELOCITY DISTRIBUTION FUNCTIONS CREATION CONTROL **************")')
         READ (9, '(A77)') buf ! --dddddd.ddd- Maximal x-velocity for e-v_x-distribution (in V_therm_e) ------")')
         READ (9, '(2x,f10.3)') Ve_x_max
         READ (9, '(A77)') buf ! --dddddd.ddd- Maximal y,z-velocity for e-v_y,z-distribution (in V_therm_e) --")')
         READ (9, '(2x,f10.3)') Ve_yz_max
         READ (9, '(A77)') buf ! ------dd----- Number of breakpoints (no evdfs if<0, all system used if 0) ---")') Creation of local distr. fun-s can be cancelled
         READ (9, '(6x,i2)') N_of_breakpoints
         READ (9, '(A77)') buf ! --dddddd.ddd- Breakpoints (ascend., node number if>0 or millimeters if <0) --")')

         DO i = 1, N_of_breakpoints
            READ (9, '(2x,f10.3)') X_breakpoint
            IF (X_breakpoint.LT.0.0_8) X_breakpoint = 0.001_8 * ABS(X_breakpoint) / delta_x_m   ! if negative - convert from millimeters
            IF (INT(X_breakpoint).EQ.0)  CYCLE                          ! must include at least one cell
            IF (INT(X_breakpoint).GE.N_cells) CYCLE                     ! must be less than the right plasma boundary to avoid ambiguity
            N_of_all_vdf_locs = N_of_all_vdf_locs + 1
            breakpoint(N_of_all_vdf_locs) = INT(X_breakpoint)
         END DO
! Save the right plasma boundary as the last breakpoint, if creation of local distribution functions was not cancelled completely
         IF (N_of_breakpoints.GE.0) THEN
            N_of_all_vdf_locs = N_of_all_vdf_locs + 1
            breakpoint(N_of_all_vdf_locs) = N_cells
         END IF

         READ (9, '(A77)') buf ! ************** ENERGY DISTRIBUTION FUNCTIONS CREATION CONTROL ***************")')
         READ (9, '(A77)') buf ! --dddddd.ddd- Maximal energy for electron e-distribution (in eV) ------------")')
         READ (9, '(2x,f10.3)') Ee_max_eV
         READ (9, '(A77)') buf ! --dddddd.ddd- Maximal energy for ion e-distribution (in eV) -----------------")')
         READ (9, '(2x,f10.3)') Ei_max_eV
         READ (9, '(A77)') buf ! ------dd----- Number of breakpoints (no edfs if<0, all system used if 0) ----")')
         READ (9, '(6x,i2)') N_of_E_breakpoints
         READ (9, '(A77)') buf ! --dddddd.ddd- Breakpoints (ascend., node number if>0 or millimeters if <0) --")')

         DO i = 1, N_of_E_breakpoints
            READ (9, '(2x,f10.3)') X_breakpoint
            IF (X_breakpoint.LT.0.0_8) X_breakpoint = 0.001_8 * ABS(X_breakpoint) / delta_x_m   ! if negative - convert from millimeters
            IF (INT(X_breakpoint).EQ.0)  CYCLE                          ! must include at least one cell
            IF (INT(X_breakpoint).GE.N_cells) CYCLE                     ! must be less than the right plasma boundary to avoid ambiguity
            N_of_all_edf_locs = N_of_all_edf_locs + 1
            breakpoint_E(N_of_all_edf_locs) = INT(X_breakpoint)
         END DO
! Save the right plasma boundary as the last breakpoint, if creation of local distribution functions was not cancelled completely
         IF (N_of_breakpoints.GE.0) THEN
            N_of_all_edf_locs = N_of_all_edf_locs + 1
            breakpoint_E(N_of_all_edf_locs) = N_cells
         END IF

         READ (9, '(A77)') buf ! --dddddd.ddd- Maximal energy for ion wall e-distribution (in eV) ------------")')
         READ (9, '(2x,f10.3)') Ei_wall_max_eV
         READ (9, '(A77)') buf ! -----ddd----- Number of bins for e-distributions (>0) -----------------------")')
         READ (9, '(5x,i3)') N_E_bins
         READ (9, '(A77)') buf ! ********************** DISTRIBUTION FUNCTION SWITCHES ***********************")')
         READ (9, '(A77)') buf ! --d--d--d--d- Velocity ( e-vx, e-vy, e-vz, i-vx, | 1/0 = on/off ) -----------")')
         READ (9, '(4(2x,i1))') flag_df_temp(1,1), flag_df_temp(1,2), flag_df_temp(1,3), flag_df_temp(1,4)
         READ (9, '(A77)') buf ! --d--d--d--d- Energy ( e, i, i l-wall, i r-wall | 1/0 = on/off ) ------------")')
         READ (9, '(4(2x,i1))') flag_df_temp(2,1), flag_df_temp(2,2), flag_df_temp(2,3), flag_df_temp(2,4)

      END IF

   ELSE

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

! report about the status of creation of probes for time dependencies ================================
! note, only the server node needs to keep the probe locations
   IF (Rank_of_process.EQ.0) THEN
      IF (N_of_probes.LE.0) THEN
         PRINT '(2x,"Probes for time dependencies are NOT specified...")'
      ELSE
         PRINT '(2x,"Probes for time dependencies are specified in ",i3," locations")', N_of_probes
! allocate the array of probe locations
         ALLOCATE(probe_node(1:N_of_probes), STAT=ALLOC_ERR)
         IF(ALLOC_ERR.NE.0)THEN
            PRINT '(/2x,"Error in ALLOCATE probe_node !!!")'
            PRINT  '(2x,"Program will be terminated now :(")'
            STOP
         END IF
         probe_node(1:N_of_probes) = probe_location(1:N_of_probes)
! write locations to the file
         OPEN (41, FILE = '_probelocs.dat')
!                        ---nnn------nnnnnn----nnnnn.nnnn
         WRITE (41, '("# number    x[node]     x[mm]")')
         DO i = 1, N_of_probes
            WRITE (41, '(3x,i3,6x,i6,4x,f10.4)') i, probe_node(i), 1000.0_8 * probe_node(i) * delta_x_m
         END DO
         CLOSE (41, STATUS = 'KEEP')
      END IF
   END IF

! report about the general status of snapshot creation
   IF (N_of_all_snaps.EQ.0) THEN

      IF (Rank_of_process.EQ.0) PRINT '(/2x,"Snapshots will NOT be created ...")'

   ELSE ! if snapshots will be created ...

! allocate the array of moments of snapshots
      IF (.NOT.ALLOCATED(Tcntr_snapshot)) THEN
         ALLOCATE(Tcntr_snapshot(1:N_of_all_snaps), STAT=ALLOC_ERR)
         IF(ALLOC_ERR.NE.0)THEN
            PRINT '(/2x,"Process ",i3," : Error in ALLOCATE Tcntr_snapshot !!!")', Rank_of_process
            PRINT  '(2x,"Program will be terminated now :(")'
            STOP
         END IF
      END IF

! move the calculated snapshot moments from the temporary array to the allocated array
      Tcntr_snapshot(1:N_of_all_snaps) = timestep(1:N_of_all_snaps)

! for server process only:
      IF (Rank_of_process.EQ.0) THEN
! save the moments to the file
         PRINT '(/2x,"The program will create ",i4," snapshots")', N_of_all_snaps

         ! IF (N_of_all_snaps.GT.0) THEN
         !    ! Print out a list of all snapshot times
         !    DO i = 1, N_of_all_snaps
         !       PRINT '(4x,"Snapshot ",i4," at time ",f10.3," ns")', i, timestep(i) * 1.0e9 * delta_t_s
         !    END DO
         ! END IF

! write moments of snapshot creation into the file
         OPEN (41, FILE = '_snapmoments.dat')
!                    "---***-----*****.*****----********---"
         WRITE (41, '("# number     time(ns)      T_cntr")')
         DO i = 1, N_of_all_snaps
            WRITE (41, '(2x,i4,5x,f11.5,4x,i8)') i, Tcntr_snapshot(i) * 1.0e9 * delta_t_s, Tcntr_snapshot(i)
         END DO
         CLOSE (41, STATUS = 'KEEP')

! report about the status of phase plane creation
         IF (Flag_pp_creation.EQ.0) THEN
            PRINT '(2x,"No phase planes will be created in snapshots ...")'
         ELSE IF (Flag_pp_creation.EQ.1) THEN
            PRINT '(2x,"Phase planes will be created: only one of each ",i5," particles will be used")', N_to_skip
         END IF

      END IF

! report about the status of local distribution functions creation ===========================
      ! First convert the temporary flags to logicals
      if (flag_df_temp(1,1).eq.0) then
         flag_evxdf = .false.
      else
         flag_evxdf = .true.
      end if
      if (flag_df_temp(1,2).eq.0) then
         flag_evydf = .false.
      else
         flag_evydf = .true.
      end if
      if (flag_df_temp(1,3).eq.0) then
         flag_evzdf = .false.
      else
         flag_evzdf = .true.
      end if
      if (flag_df_temp(1,4).eq.0) then
         flag_ivxdf = .false.
      else
         flag_ivxdf = .true.
      end if

      ! Next report the status to the output file
      if (flag_evxdf.or.flag_evydf.or.flag_evzdf.or.flag_ivxdf) then 
         IF (N_of_all_vdf_locs.EQ.0) THEN
            IF (Rank_of_process.EQ.0) PRINT '(2x,"No local velocity distribution functions will be created in snapshots ...")'
         ELSE
            IF (Rank_of_process.EQ.0) PRINT '(2x,"The velocity distribution functions will be calculated in ",i3," locations")', &
            & N_of_all_vdf_locs
            ! allocate the array of boundaries of locations for calculation of VDF in snapshots
            IF (.NOT.ALLOCATED(Vdf_location_bnd)) THEN
               ALLOCATE(Vdf_location_bnd(1:N_of_all_vdf_locs), STAT=ALLOC_ERR)
               IF(ALLOC_ERR.NE.0)THEN
                  PRINT '(/2x,"Process ",i3," : Error in ALLOCATE Vdf_location_bnd !!!")', Rank_of_process
                  PRINT  '(2x,"Program will be terminated now :(")'
                  STOP
               END IF
            END IF
            ! move the calculated locations for VDF from the temporary array to the allocated array
            Vdf_location_bnd(1:N_of_all_vdf_locs) = breakpoint(1:N_of_all_vdf_locs)
            ! for server process only:
            IF (Rank_of_process.EQ.0) THEN
               ! write locations to the file
               OPEN (41, FILE = '_vdflocations.dat')
               !            "---***---#*.*****E#**---#*.*****E#**------******------******"
               WRITE (41, '("# number   left(mm)       right(mm)     left node   right node")')
               WRITE (41, '(3x,i3,2(3x,e12.5),2(6x,i6))') 1, &
               & 0.0, &
               & 1000.0_8 * Vdf_location_bnd(1) * delta_x_m, &
               & 0, &
               & Vdf_location_bnd(1)
               DO i = 2, N_of_all_vdf_locs
                  WRITE (41, '(3x,i3,2(3x,e12.5),2(6x,i6))') i, &
                  & 1000.0_8 * Vdf_location_bnd(i-1) * delta_x_m, &
                  & 1000.0_8 * Vdf_location_bnd(i) * delta_x_m, &
                  & Vdf_location_bnd(i-1), &
                  & Vdf_location_bnd(i)
               END DO
               CLOSE (41, STATUS = 'KEEP')
            END IF
         END IF
      end if

      N_box_Vx_e      = Ve_x_max * N_box_vel - 1
      N_box_Vyz_e     = Ve_yz_max * N_box_vel - 1
      N_box_Vx_e_low  = -N_box_Vx_e                   !
      N_box_Vx_e_top  =  N_box_Vx_e + 1               ! this will save three operations of addition and
      N_box_Vyz_e_low = -N_box_Vyz_e                  ! three operations of sign changing for each electron
      N_box_Vyz_e_top =  N_box_Vyz_e + 1              !

      IF (N_box_Vyz_e.GT.N_box_Vx_e) THEN
         N_box_Vx_i = N_box_Vyz_e
      ELSE
         N_box_Vx_i = N_box_Vx_e
      END IF
      N_box_Vx_i_low  = -N_box_Vx_i                   ! this will save one addition and one sign changing for each ion
      N_box_Vx_i_top  =  N_box_Vx_i + 1               !

! report about the status of local energy distribution functions creation ===========================
      ! First convert the temporary flags to logicals
      if (flag_df_temp(2,1).eq.0) then
         flag_eedf = .false.
      else
         flag_eedf = .true.
      end if
      if (flag_df_temp(2,2).eq.0) then
         flag_iedf = .false.
      else
         flag_iedf = .true.
      end if
      if (flag_df_temp(2,3).eq.0) then
         flag_ilwedf = .false.
      else
         flag_ilwedf = .true.
      end if
      if (flag_df_temp(2,4).eq.0) then
         flag_irwedf = .false.
      else
         flag_irwedf = .true.
      end if

      ! Next report the status to the output file
      if (flag_eedf.or.flag_iedf) then 
         IF (N_of_all_edf_locs.EQ.0) THEN
            IF (Rank_of_process.EQ.0) PRINT '(2x,"No local energy distribution functions will be created in snapshots ...")'
         ELSE
            IF (Rank_of_process.EQ.0) PRINT '(2x,"The energy distribution functions will be calculated in ",i3," locations")', &
            & N_of_all_edf_locs
            ! allocate the array of boundaries of locations for calculation of VDF in snapshots
            IF (.NOT.ALLOCATED(Edf_location_bnd)) THEN
               ALLOCATE(Edf_location_bnd(1:N_of_all_edf_locs), STAT=ALLOC_ERR)
               IF(ALLOC_ERR.NE.0)THEN
                  PRINT '(/2x,"Process ",i3," : Error in ALLOCATE Edf_location_bnd !!!")', Rank_of_process
                  PRINT  '(2x,"Program will be terminated now :(")'
                  STOP
               END IF
            END IF
            ! move the calculated locations for EDF from the temporary array to the allocated array
            Edf_location_bnd(1:N_of_all_edf_locs) = breakpoint_E(1:N_of_all_edf_locs)
            ! for server process only:
            IF (Rank_of_process.EQ.0) THEN
               ! write locations to the file
               OPEN (41, FILE = '_edflocations.dat')
               !            "---***---#*.*****E#**---#*.*****E#**------******------******"
               WRITE (41, '("# number   left(mm)       right(mm)     left node   right node")')
               WRITE (41, '(3x,i3,2(3x,e12.5),2(6x,i6))') 1, &
               & 0.0, &
               & 1000.0_8 * Edf_location_bnd(1) * delta_x_m, &
               & 0, &
               & Edf_location_bnd(1)
               DO i = 2, N_of_all_edf_locs
                  WRITE (41, '(3x,i3,2(3x,e12.5),2(6x,i6))') i, &
                  & 1000.0_8 * Edf_location_bnd(i-1) * delta_x_m, &
                  & 1000.0_8 * Edf_location_bnd(i) * delta_x_m, &
                  & Edf_location_bnd(i-1), &
                  & Edf_location_bnd(i)
               END DO
               CLOSE (41, STATUS = 'KEEP')
            END IF
         END IF
      end if

      delta_Ee_eV = Ee_max_eV / N_E_bins
      delta_Ei_eV = Ei_max_eV / N_E_bins
      delta_Ei_wall_eV = Ei_wall_max_eV / N_E_bins

! allocate arrays, related with the distribution functions, initiate and write to the files arrays of middles of velocity boxes
      CALL CREATE_DF_ARRAYS

! create standard distribution functions, which can be used for benchmarking of distributions produced in simulations
      IF ((Rank_of_process.EQ.0).AND.(N_of_all_vdf_locs.GT.0)) CALL CREATE_STANDARD_DISTRIBUTIONS

   END IF

! if particles are initialized NOT from the checkpoint datafiles                  (Restore_from_checkpoint = 0)
! or if we initialize particles from the checkpoint but start from the beginning  (Restore_from_checkpoint = 2)
   IF (Restore_from_checkpoint.NE.1) THEN
      current_snap = 1
      Accumulate_wall_df = .FALSE.
   END IF

! always set this counter to zero unless in future you add reopening of device #90 and file _NNNN_evolF_vsx.dat
   counter_of_profiles_to_save = 0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

   ! Number of timesteps per RF period
   steps_per_period = 1. / (f_rf_Hz * delta_t_s)

   ! if flag_Out is negative or zero then WriteOut_step is given in units of RF periods
   if (flag_Out.le.0) then
      if (.not.rf_on) then 
         print '(/2x,"Process ",i3," : WriteOut_Step in ssc_diagnostics.dat was given in terms of RF period")', Rank_of_process
         print '(/16x,"but RF is turned off. Change method of specifying WriteOut_Step to continue.")'
         print '(2x,"Program will be terminated now")'
         stop
      endif
      if (Rank_of_process.EQ.0) print '(/2x,"Data output interval was given in units of RF period ...")'
      ! Multiply the output step (expanded by the flag) by the number of steps per RF period
      NumOut_RFper = WriteOut_step * (10.**flag_Out)
      WriteOut_step = NumOut_RFper * steps_per_period
   else
      ! if WriteOut_step is negative then it is given in plasma periods, assumed to be integer !!!
      IF (WriteOut_step.LT.0) THEN
         IF (Rank_of_process.EQ.0) PRINT '(/2x,"Data output interval was given in units of electron plasma period ...")'
         WriteOut_step = ABS(WriteOut_step) * 6.28318530718_8 / (W_plasma_s1 * delta_t_s) 
      END IF
   endif

   ! if flag_Start is negative or zero then WriteStart_step is given in units of RF periods
   if (flag_Start.le.0) then
      if (.not.rf_on) then 
         print '(/2x,"Process ",i3," : WriteStart_Step in ssc_diagnostics.dat was given in terms of RF period")', Rank_of_process
         print '(/16x,"but RF is turned off. Change method of specifying WriteStart_Step to continue.")'
         print '(2x,"Program will be terminated now")'
         stop
      endif
      if (Rank_of_process.EQ.0) print '(/2x,"Data collection start moment was given in units of RF period ...")'
      ! Multiply the start step (expanded by the flag) by the number of steps per RF period
      NumStart_RFper = WriteStart_step * (10.**flag_Start)
      WriteStart_step = NumStart_RFper * steps_per_period
   else
      ! if WriteStart_step is negative then it is given in plasma periods, assumed to be integer !!!
      IF (WriteStart_step.LT.0) THEN
         IF (Rank_of_process.EQ.0) PRINT '(/2x,"Data collection start moment was given in units of electron plasma period ...")'
         WriteStart_step = ABS(WriteStart_step) * 6.28318530718_8 / (W_plasma_s1 * delta_t_s)
      END IF
   endif

   ! if flag_Avg is negative or zero then WriteAvg_step is given in units of RF periods
   if (flag_Avg.le.0) then
      if (.not.rf_on) then 
         print '(/2x,"Process ",i3," : WriteAvg_Step in ssc_diagnostics.dat was given in terms of RF period")', Rank_of_process
         print '(/16x,"but RF is turned off. Change method of specifying WriteAvg_Step to continue.")'
         print '(2x,"Program will be terminated now")'
         stop
      endif
      if (Rank_of_process.EQ.0) print '(2x,"Data averaging window was given in units of RF period ...")'
      ! Multiply the averaging step (expanded by the flag) by the number of steps per RF period
      NumAvg_RFper = WriteAvg_step * (10.**flag_Avg)
      WriteAvg_step = NumAvg_RFper * steps_per_period
   else
      ! if WriteAvg_step is negative then it is given in plasma periods, assumed to be integer !!!
      IF (WriteAvg_step.LT.0) THEN
         IF (Rank_of_process.EQ.0) PRINT '(2x,"Data averaging window was given in units of electron plasma period ...")'
         WriteAvg_step = ABS(WriteAvg_step) * 6.28318530718_8 / (W_plasma_s1 * delta_t_s)
      END IF
   endif

   if (WriteAvg_step.gt.WriteOut_step) then
      if (Rank_of_process.EQ.0) then
         print '(/2x,"Process ",i3," : Data averaging window (WriteAvg_step) is larger than data output interval (WriteOut_step) !!!")', Rank_of_process
         print  '(2x,"Program will be terminated now :(")'
      endif
      stop
   endif

! ensure that SaveCheck_step [defined previously in INITIATE_PARAMETERS] is a multiple of WriteOut_step
!** but never zero in case checkpoints were requested ** Oct. 2016
   if (SaveCheck_step.gt.0) then
      SaveCheck_step = (SaveCheck_step / WriteOut_Step) * WriteOut_step
      if (SaveCheck_step.eq.0) SaveCheck_step = WriteOut_step
   else
   endif

   IF (Rank_of_process.EQ.0) THEN
      PRINT '(/2x,"Diagnostic values will be written into the file with time interval : ",f13.4," ns")', WriteOut_step * delta_t_s * 1.0e9
      if (flag_Out.le.0) print '(67x,"or, ",f9.4," RF periods")', NumOut_RFper
      PRINT  '(2x,"First data will start collecting at the moment                     : ",f13.4," ns")',  WriteStart_step * delta_t_s * 1.0e9
      if (flag_Start.le.0) print '(67x,"or, ",f9.4," RF periods")', NumStart_RFper
      PRINT  '(2x,"Averaging will be carried out during the following time interval   : ",f13.4," ns")',   WriteAvg_step * delta_t_s * 1.0e9
      if (flag_Start.le.0) print '(67x,"or, ",f9.4," RF periods")', NumAvg_RFper
      IF (SaveCheck_step.GT.0) THEN
         PRINT '(2x,"Checkpoints will be created with time interval                    : ",f13.4," ns")',  SaveCheck_step * delta_t_s * 1.0e9
      ELSE
         PRINT '(2x,"Checkpoints will NOT be created...")'
      END IF
   END IF

   CALL AllocateDiagnosticArrays

! if particles are initialized NOT from the checkpoint datafiles                  (Restore_from_checkpoint = 0)
! or if we initialize particles from the checkpoint but start from the beginning  (Restore_from_checkpoint = 2)
   IF (Restore_from_checkpoint.NE.1) Start_diag_Tcntr = WriteStart_step             ! ### ###
   Finish_diag_Tcntr = Start_diag_Tcntr + WriteAvg_step - 1

! clear the diagnostics arrays
   VX2_cell = 0.0_8
   VY2_cell = 0.0_8
   VZ2_cell = 0.0_8

   VX_cell = 0.0_8
   VY_cell = 0.0_8
   VZ_cell = 0.0_8

   Npart_cell = 0

   QF_mesh = 0.0_8

   QVX_mesh = 0.0_8
   QVY_mesh = 0.0_8
   QVZ_mesh = 0.0_8

   N_new_cell = 0 !integer
   P_heat_cell = 0.0_8
   NVX_mesh = 0.0_8

! clear the dim-less quantities accumulated only during the last timestep of the averaging period
   Energy_pot       = 0.0_8
   Avg_kin_energy_x = 0.0_8           ! array
   Avg_kin_energy_y = 0.0_8           ! array
   Avg_kin_energy_z = 0.0_8           ! array
   VY_recent        = 0.0_8           ! array
   VZ_recent        = 0.0_8           ! array

! clear the energy rates for collisions and heating (will be accumulated)
   Rate_energy_coll = 0.0_8           ! array
   Rate_energy_heat = 0.0_8

! clear the numbers/energies accumulated on different walls in (final/pre push, emission)
   Rate_number_leftwall  = 0          ! array
   Rate_number_rightwall = 0          ! array
   Rate_number_leftemit  = 0          ! array
   Rate_number_rightemit = 0          ! array
   Rate_energy_leftwall  = 0.0_8      ! array
   Rate_energy_rightwall = 0.0_8      ! array
   Rate_energy_leftemit  = 0.0_8      ! array
   Rate_energy_rightemit = 0.0_8      ! array

! calculate the factors
   Averaging_factor        = 1.0_8 / WriteAvg_step
   J_scl_Am2               = Averaging_factor * (N_plasma_m3 / DBLE(N_of_particles_cell)) * e_Cl * (v_Te_ms / DBLE(N_box_vel))
   N_in_macro              = N_plasma_m3 * delta_x_m / dble(N_of_particles_cell)
   Factor_energy_eV        = T_e_eV / N_box_vel**2
   Factor_energy_macro_eV  = N_in_macro * Factor_energy_eV
   Factor_energy_pot_eV    = N_in_macro * (2.0_8 * T_e_eV / r_cells_debye) * 0.5_8
   if (Factor_energy_pot_eV .eq. 0.) write(*,*) 'XOXOXO'
   Factor_rate_ns1         = 1.0_8 / (WriteOut_step * delta_t_s * 1.0e9)
   Factor_rate_macro_ns1   = N_in_macro * Factor_rate_ns1
   Factor_rate_eVns1       = Factor_energy_eV * Factor_rate_macro_ns1
   Factor_Joule_energy_eV  = N_in_macro * (v_Te_ms / N_box_vel) * E_z_ext_Vm * (delta_t_s)

   f_factor = 0.0_8      ! array (1:N_spec) set zero here in order to avoid problems when N_spec = 1 (i.e. f_factor(2)=0.0_8 always in this case)

!  Factor_energy_eV  = 2.0_8 * T_e_eV * N_plasma_m3 * L_debye_m / N_of_particles_debye
!  Factor_rate_eVns1 = Factor_energy_eV / (WriteAvg_step * delta_t_s * 1.0e9)

   Save_check_Tcntr = Finish_diag_Tcntr + SaveCheck_step

   print '("Process ",i4," : finished INITIATE_DIAGNOSTICS")', Rank_of_process

END SUBROUTINE INITIATE_DIAGNOSTICS

!------------------------------
subroutine CALCULATE_SNAPSHOT_TIMES(Aprx_snap_start_ns, Aprx_snap_finish_ns, Aprx_n_of_snaps, T2_old, N2_old, timestep)
!! Calculates the exact snapshot output timesteps for each snapshot group requested
!! in ssc_diagnostics.dat. Each snapshot timestep coincides with a diagnostic
!! Finish_Diag_Tcntr timestep.
!!
!! Because Start_diag_Tcntr steps (and thus Finish_diag_Tcntr steps) are spaced apart
!! by jumps of length WriteOut_step, this subroutine calculates the first snapshot timestep
!! as a multiple of WriteOut_step after the first diagnostic output timestep and places
!! all remaining snapshots within the group evenly spaced between the first timestep and 
!! the final target timestep (Aprx_snap_finish_ns).
!!
!! Outside the program we read the parameters of current set of snapshot from the data file:
!!
!!   READ (9, '(12x,f10.3,12x,f10.3,14x,i4)') Aprx_snap_start_ns, Aprx_snap_finish_ns, Aprx_n_of_snaps

   use CurrentProblemValues, only : delta_t_s, Max_T_cntr
   use Diagnostics, only : WriteOut_step, WriteStart_step, WriteAvg_step
   use Snapshots, only : N_of_all_snaps

   IMPLICIT NONE
   real(8), intent(in) :: Aprx_snap_start_ns       !! Requested time [ns] of first snapshot output, passed from file
   real(8), intent(in) :: Aprx_snap_finish_ns      !! Requested time [ns] of final snapshot output, passed from file
   integer, intent(in) :: Aprx_n_of_snaps          !! Requested number of snapshots in group, passed from file
   integer, intent(inout) :: timestep(1:9999)      !! Array of snapshot timesteps
   integer, intent(inout) :: T2_old                !! Finish timestep for prior snapshot group
   integer, intent(inout) :: N2_old                !! Number of WriteOut_steps between the first diagnostic output and T2_old

   integer T1                                      !! Start timestep for snapshot group
   integer T2                                      !! End/output timestep for snapshot group
   integer N1                                      !! Number of WriteOut_steps between the first diagnostic output and T1
   integer N2                                      !! Number of WriteOut_steps between the first diagnostic output and T2

   integer first_diag_output_step                  !! Temporary variable to increase readability of this function.
                                                   !!   This is calculated later as a module variable and named Finish_diag_Tcntr
   integer n                                       !! Snapshot index
   integer Fact_n_of_snaps                         !! Actual number of snapshots in group
   integer large_step                              !! Number of WriteOut_steps between T1 and T2 dived by the requested number
                                                   !!   of snapshots in the group. If large_step >=1, then the requested number
                                                   !!   is the actual number of snapshots in the group.
                                                   !! This factor helps ensure snapshots are evenly spaced within the group. If
                                                   !!   1 <= large_step < 2 when caculated, then snapshots are outputted every
                                                   !!  WriteOut_step. If 2 <= large_step < 3, then snapshots are outputted every
                                                   !!  second WriteOut_step, etc.

   ! Skip if zero snapshots were requested in the group
   IF (Aprx_n_of_snaps.LT.1) RETURN

   first_diag_output_step = WriteStart_step + WriteAvg_step - 1

   ! Calculate group start timestep asked for in the file
   T1 = Aprx_snap_start_ns / (delta_t_s * 1.0e9)
   IF (T1.LT.first_diag_output_step) T1 = first_diag_output_step ! If T1 < first_diag_output_step, set T1 = first_diag_output_step
   N1 = (T1 - first_diag_output_step) / WriteOut_step            ! Calculate number of write steps between T1 and first_diag_output_step
   T1 = first_diag_output_step + N1 * WriteOut_step              ! Make T1 coincide with a write step after first_diag_output_step

   ! Calculate group end timestep asked for in the file
   T2 = Aprx_snap_finish_ns / (delta_t_s * 1.0e9)
   IF (T2.LT.first_diag_output_step) T2 = first_diag_output_step
   N2 = (T2 - first_diag_output_step) / WriteOut_step
   T2 = first_diag_output_step + N2 * WriteOut_step              ! Make T2 coincide with a write step after first_diag_output_step

   ! Ensure T1 falls after T2 for the previous group of snapshots
   IF (T1.LE.T2_old) THEN
      T1 = T2_old + WriteOut_step
      N1 = N2_old + 1
   END IF

   ! Ensure the T2 falls before the last simulation timestep
   DO WHILE (T2.GT.Max_T_cntr)
      T2 = T2 - WriteOut_step
      N2 = N2 - 1
   END DO

   ! Skip this group if T1 falls after T2
   IF (T1.GT.T2) RETURN

   ! ****************
   ! If we are here than the T1 and T2 can be used for calculation of moments of snapshots
   ! ****************

   ! Calculate the actual number of snapshots in the current group:
   IF (Aprx_n_of_snaps.EQ.1) THEN ! If only one snapshot is requested:
      large_step = 0
      Fact_n_of_snaps = 1
      T2 = T1  ! make the snapshot output at T1
      N2 = N1

   ELSE ! If more than one snapshot, calculate the number of WriteOut_steps between the first and last snapshots
      ! Determine the number of WriteOut_steps between T1 and T2 (if large_step >=1, then the requested number of snapshots can be created)
      large_step = (N2 - N1) / (Aprx_n_of_snaps - 1)

      IF (large_step.EQ.0) THEN
         large_step = 1
         Fact_n_of_snaps = N2 - N1 + 1 ! If the number of snapshots requested is larger than the number of
                                       !   WriteOut_steps available, then make the number of snapshots equal
                                       !   to the number of WriteOut_steps between T1 and T2 (including an 
                                       !   initial snapshot output at T1)
      ELSE
         Fact_n_of_snaps = Aprx_n_of_snaps                ! There are over enough WriteOut_steps to make the requested number of snapshots
         N2 = N1 + large_step * (Fact_n_of_snaps - 1)     ! Calculate the number of WriteOut_steps between first_diag_output_step
                                                          !   and the last snapshot of the group
         T2 = first_diag_output_step + N2 * WriteOut_step ! Calculate the last snapshot output timestep
      END IF
   END IF

   ! for each (possible) snapshot in the group save the snapshot output moment in the temporary array
   DO n = 1, Fact_n_of_snaps
      N_of_all_snaps = N_of_all_snaps + 1
      timestep(N_of_all_snaps) =  T1 + (n - 1) * large_step * WriteOut_step
   END DO

   ! save the finish data for the current set of snapshots
   T2_old = T2
   N2_old = N2

end subroutine CALCULATE_SNAPSHOT_TIMES

!------------------------------
!
SUBROUTINE PREPARE_TIME_DEPENDENCE_DATAFILES

   USE ParallelOperationValues
   USE CurrentProblemValues
   USE Diagnostics
   IMPLICIT NONE

   LOGICAL exists
   INTEGER i
!  INTEGER i_dummy
   REAL    r_dummy, arr_dummy(1:99)

   INTEGER temp_N_of_records

   IF (Rank_of_process.NE.0) RETURN

! if particles are initialized NOT from the checkpoint datafiles                  (Restore_from_checkpoint = 0)
! or if we initialize particles from the checkpoint but start from the beginning  (Restore_from_checkpoint = 2)
   IF (Restore_from_checkpoint.NE.1) THEN
! to distinguish between the first calculation of energies and the consecutive calculations:
      Energy_full_eV = 0.0_8
      Energy_pot_eV  = 0.0_8
      Energy_kin_eV  = 0.0_8             ! array
      Energy_wall_eV = 0.0_8             ! array
      Energy_emit_eV = 0.0_8             ! array
      Energy_coll_eV = 0.0_8             ! array
      Energy_heat_eV = 0.0_8
      Enr_consumed_Jm2 = 0.0_8           !updated only for capacitive discharge
      Sigma_electrode = 0.0_8            !same

! initiate counter for text output skipping. Text is printed when the very first time diagnostics is obtained
! and then each TextOut_avg_prds_to_skip diagnostics are performed without text output.
      text_output_counter = TextOut_avg_prds_to_skip
! initiate the number of records saved in the data file (temporal evolution)
      N_of_saved_records = 0

! open/close files in order to replace previously existed files
! ACCUMULATED ENERGY VALUES

      OPEN (50, FILE = 'dim_fullenergy_vst.dat', STATUS = 'REPLACE')          ! full energy of system (kinetic + potential)
      CLOSE (50, STATUS = 'KEEP')

      OPEN (50, FILE = 'dim_potenergy_vst.dat', STATUS = 'REPLACE')           ! potential energy of system
      CLOSE (50, STATUS = 'KEEP')

      OPEN (50, FILE = 'dim_kinenergy_vst.dat', STATUS = 'REPLACE')           ! kinetic energy of system: all species / electrons / ions
      CLOSE (50, STATUS = 'KEEP')

      OPEN (50, FILE = 'dim_energy_1e_vst.dat', STATUS = 'REPLACE')           ! average kinetic energy of one electron: V / X / Y / Z
      CLOSE (50, STATUS = 'KEEP')

      OPEN (50, FILE = 'dim_energy_1i_vst.dat', STATUS = 'REPLACE')           ! average kinetic energy of one ion:      V / X / Y / Z
      CLOSE (50, STATUS = 'KEEP')

      OPEN (50, FILE = 'dim_energy_wall_vst.dat', STATUS = 'REPLACE')         ! kinetic energy of particles that have struck the walls:
      CLOSE (50, STATUS = 'KEEP')                                             ! all species / electrons / ions

      OPEN (50, FILE = 'dim_energy_emit_vst.dat', STATUS = 'REPLACE')         ! kinetic energy of particles that were emitted/reflected from the walls:
      CLOSE (50, STATUS = 'KEEP')                                             ! all species / electrons / ions

      OPEN (50, FILE = 'dim_energy_coll_vst.dat', STATUS = 'REPLACE')         ! kinetic energy, lost due to collisions with neutrals:
      CLOSE (50, STATUS = 'KEEP')                                             ! all species / electrons / ions

      OPEN (50, FILE = 'dim_energy_heat_vst.dat', STATUS = 'REPLACE')         ! energy of Joule heating of plasma by accelerating electric field
      CLOSE (50, STATUS = 'KEEP')

      OPEN (50, FILE = 'dim_energy_cons_violat_vst.dat', STATUS = 'REPLACE')  ! violation of conservation of energy (%)
      CLOSE (50, STATUS = 'KEEP')

! RATES OF CHANGE OF ENERGY VALUES

      OPEN (50, FILE = 'dim_rate_fullenergy_vst.dat', STATUS = 'REPLACE')     ! full energy of system (kinetic + potential)
      CLOSE (50, STATUS = 'KEEP')

      OPEN (50, FILE = 'dim_rate_potenergy_vst.dat', STATUS = 'REPLACE')      ! potential energy of system
      CLOSE (50, STATUS = 'KEEP')

      OPEN (50, FILE = 'dim_rate_kinenergy_vst.dat', STATUS = 'REPLACE')      ! kinetic energy of system: all species / electrons / ions
      CLOSE (50, STATUS = 'KEEP')

      OPEN (50, FILE = 'dim_rate_energy_wall_vst.dat', STATUS = 'REPLACE')    ! kinetic energy of particles that have struck the walls:
      CLOSE (50, STATUS = 'KEEP')                                             ! all species / electrons / ions

      OPEN (50, FILE = 'dim_rate_energy_emit_vst.dat', STATUS = 'REPLACE')    ! kinetic energy of particles that were emitted/reflected from the walls:
      CLOSE (50, STATUS = 'KEEP')                                             ! all species / electrons / ions

      OPEN (50, FILE = 'dim_rate_energy_coll_vst.dat', STATUS = 'REPLACE')    ! kinetic energy, lost due to collisions with neutrals:
      CLOSE (50, STATUS = 'KEEP')                                             ! all species / electrons / ions

      OPEN (50, FILE = 'dim_rate_energy_heat_vst.dat', STATUS = 'REPLACE')    ! energy of Joule heating of plasma by accelerating electric field
      CLOSE (50, STATUS = 'KEEP')

! FLOW VALUES, AVERAGED OVER CROSS SECTION

      OPEN (50, FILE = 'dim_vy_e_vst.dat', STATUS = 'REPLACE')                ! velocity of electron Y-flow
      CLOSE (50, STATUS = 'KEEP')

      OPEN (50, FILE = 'dim_vz_e_vst.dat', STATUS = 'REPLACE')                ! velocity of electron Z-flow
      CLOSE (50, STATUS = 'KEEP')

      OPEN (50, FILE = 'dim_jy_vst.dat', STATUS = 'REPLACE')                  ! electric current density in Y-direction
      CLOSE (50, STATUS = 'KEEP')

      OPEN (50, FILE = 'dim_jz_vst.dat', STATUS = 'REPLACE')                  ! electric current density in Z-direction
      CLOSE (50, STATUS = 'KEEP')

! PLASMA DENSITIES, POTENTIAL, ELECTRIC FIELD, IN PROBES

      OPEN (50, FILE = 'dim_Ne_vst.dat', STATUS = 'REPLACE')                  ! electrons
      CLOSE (50, STATUS = 'KEEP')

      OPEN (50, FILE = 'dim_Ni_vst.dat', STATUS = 'REPLACE')                  ! ions
      CLOSE (50, STATUS = 'KEEP')

      OPEN (50, FILE = 'dim_F_vst.dat', STATUS = 'REPLACE')                   ! middle
      CLOSE (50, STATUS = 'KEEP')

      OPEN (50, FILE = 'dim_Ex_vst.dat', STATUS = 'REPLACE')                  ! wall: left / right
      CLOSE (50, STATUS = 'KEEP')

! COLLISION FREQUENCIES

      OPEN (50, FILE = 'dim_fen_collisions_vst.dat', STATUS = 'REPLACE')      ! frequencies of electron-neutral collisions: elastic / excitation-1 / excitation-2 / ionization / turbulence
      CLOSE (50, STATUS = 'KEEP')

      OPEN (50, FILE = 'dim_fin_collisions_vst.dat', STATUS = 'REPLACE')      ! frequencies of ion-neutral collisions: elastic / charge exchange / turbulence
      CLOSE (50, STATUS = 'KEEP')

! AT THE WALLS

      OPEN (50, FILE = 'dim_n_leftwall_vst.dat', STATUS = 'REPLACE')          ! number of particles that hit the left wall: electrons / ions
      CLOSE (50, STATUS = 'KEEP')

      OPEN (50, FILE = 'dim_n_rightwall_vst.dat', STATUS = 'REPLACE')         ! number of particles that hit the right wall: electrons / ions
      CLOSE (50, STATUS = 'KEEP')

      OPEN (50, FILE = 'dim_avgenergy_leftwall_vst.dat', STATUS = 'REPLACE')  ! average energy of particles that hit the left wall: electrons / ions
      CLOSE (50, STATUS = 'KEEP')

      OPEN (50, FILE = 'dim_avgenergy_rightwall_vst.dat', STATUS = 'REPLACE') ! average energy of particles that hit the right wall: electrons / ions
      CLOSE (50, STATUS = 'KEEP')

      OPEN (50, FILE = 'dim_avgenergy_leftemit_vst.dat', STATUS = 'REPLACE')  ! average energy of particles emitted from the left wall: electrons / ions
      CLOSE (50, STATUS = 'KEEP')

      OPEN (50, FILE = 'dim_avgenergy_rightemit_vst.dat', STATUS = 'REPLACE') ! average energy of particles emitted from the right wall: electrons / ions
      CLOSE (50, STATUS = 'KEEP')

      OPEN (50, FILE = 'dim_gamma_left_vst.dat', STATUS = 'REPLACE')          ! emission coefficient at the left wall: electrons / ions
      CLOSE (50, STATUS = 'KEEP')

      OPEN (50, FILE = 'dim_gamma_right_vst.dat', STATUS = 'REPLACE')         ! emission coefficient at the right wall: electrons / ions
      CLOSE (50, STATUS = 'KEEP')

      OPEN (50, FILE = 'dim_n_leftwall_comp_vst.dat', STATUS = 'REPLACE')          ! number of electrons that hit the left wall,
      CLOSE (50, STATUS = 'KEEP')                                                  ! components: secondary / collided / unrecognized

      OPEN (50, FILE = 'dim_avgenergy_leftwall_comp_vst.dat', STATUS = 'REPLACE')  ! average energy of el-ns that hit the left wall,
      CLOSE (50, STATUS = 'KEEP')                                                  ! components: secondary / collided / unrecognized

      OPEN (50, FILE = 'dim_gamma_left_comp_vst.dat', STATUS = 'REPLACE')          ! emission coefficient at the left wall,
      CLOSE (50, STATUS = 'KEEP')                                                  ! components: secondary / collided / unrecognized

      OPEN (50, FILE = 'dim_n_rightwall_comp_vst.dat', STATUS = 'REPLACE')         ! number of electrons that hit the right wall,
      CLOSE (50, STATUS = 'KEEP')                                                  ! components: secondary / collided / unrecognized

      OPEN (50, FILE = 'dim_avgenergy_rightwall_comp_vst.dat', STATUS = 'REPLACE') ! average energy of el-ns that hit the right wall,
      CLOSE (50, STATUS = 'KEEP')                                                  ! components: secondary / collided / unrecognized

      OPEN (50, FILE = 'dim_gamma_right_comp_vst.dat', STATUS = 'REPLACE')         ! emission coefficient at the right wall,
      CLOSE (50, STATUS = 'KEEP')                                                  ! components: secondary / collided / unrecognized

   ELSE         ! if the particles are initialized from the checkpoint datafiles     (Restore_from_checkpoint = 1)
      ! if the protocol files exist we set the length of existing protocol files equal to "N_of_saved_records" records
      ! later we will open the protocol files with " position = 'append' "
      ! therefore, if the protocol files does not exist, they will be written from the beginning

! full energy of system (kinetic + potential)
      INQUIRE (FILE = 'dim_fullenergy_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_fullenergy_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,2x,e12.5)') r_dummy !, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! potential energy of system
      INQUIRE (FILE = 'dim_potenergy_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_potenergy_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,2x,e12.5)') r_dummy !, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! kinetic energy of system: all species / electrons / ions
      INQUIRE (FILE = 'dim_kinenergy_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_kinenergy_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,3(2x,e12.5))') r_dummy !, r_dummy, r_dummy, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! average kinetic energy of one electron: V / X / Y / Z
      INQUIRE (FILE = 'dim_energy_1e_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_energy_1e_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,4(2x,e12.5))') r_dummy !, r_dummy, r_dummy, r_dummy, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! average kinetic energy of one ion:      V / X / Y / Z
      INQUIRE (FILE = 'dim_energy_1i_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_energy_1i_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,4(2x,e12.5))') r_dummy !, r_dummy, r_dummy, r_dummy, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! kinetic energy of particles that have struck the walls: all species / electrons / ions
      INQUIRE (FILE = 'dim_energy_wall_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_energy_wall_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,3(2x,e12.5))') r_dummy !, r_dummy, r_dummy, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! kinetic energy of particles that were emitted/reflected from the walls: all species / electrons / ions
      INQUIRE (FILE = 'dim_energy_emit_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_energy_emit_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,3(2x,e12.5))') r_dummy !, r_dummy, r_dummy, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! kinetic energy, lost due to collisions with neutrals: all species / electrons / ions
      INQUIRE (FILE = 'dim_energy_coll_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_energy_coll_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,3(2x,e12.5))') r_dummy !, r_dummy, r_dummy, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! energy of Joule heating of plasma by accelerating electric field
      INQUIRE (FILE = 'dim_energy_heat_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_energy_heat_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,2x,e12.5)') r_dummy !, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! violation of conservation of energy (%)
      INQUIRE (FILE = 'dim_energy_cons_violat_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_energy_cons_violat_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,2x,e12.5)') r_dummy !, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! RATES OF CHANGE OF ENERGY VALUES

! full energy of system (kinetic + potential)
      INQUIRE (FILE = 'dim_rate_fullenergy_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_rate_fullenergy_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,2x,e12.5)') r_dummy !, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! potential energy of system
      INQUIRE (FILE = 'dim_rate_potenergy_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_rate_potenergy_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,2x,e12.5)') r_dummy !, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! kinetic energy of system: all species / electrons / ions
      INQUIRE (FILE = 'dim_rate_kinenergy_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_rate_kinenergy_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,3(2x,e12.5))') r_dummy !, r_dummy, r_dummy, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! kinetic energy of particles that have struck the walls: all species / electrons / ions
      INQUIRE (FILE = 'dim_rate_energy_wall_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_rate_energy_wall_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,3(2x,e12.5))') r_dummy !, r_dummy, r_dummy, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! kinetic energy of particles that were emitted/reflected from the walls: all species / electrons / ions
      INQUIRE (FILE = 'dim_rate_energy_emit_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_rate_energy_emit_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,3(2x,e12.5))') r_dummy !, r_dummy, r_dummy, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! kinetic energy, lost due to collisions with neutrals: all species / electrons / ions
      INQUIRE (FILE = 'dim_rate_energy_coll_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_rate_energy_coll_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,3(2x,e12.5))')  r_dummy !, r_dummy, r_dummy, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! energy of Joule heating of plasma by accelerating electric field
      INQUIRE (FILE = 'dim_rate_energy_heat_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_rate_energy_heat_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,2x,e12.5)') r_dummy !, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! FLOW VALUES, AVERAGED OVER CROSS SECTION

! velocity of electron Y-flow
      INQUIRE (FILE = 'dim_vy_e_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_vy_e_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,2x,e12.5)') r_dummy !, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! velocity of electron Z-flow
      INQUIRE (FILE = 'dim_vz_e_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_vz_e_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,2x,e12.5)') r_dummy !, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! electric current density in Y-direction
      INQUIRE (FILE = 'dim_jy_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_jy_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,2x,e12.5)') r_dummy !, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! electric current density in Z-direction
      INQUIRE (FILE = 'dim_jz_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_jz_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,2x,e12.5)') r_dummy !, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! PLASMA DENSITIES, POTENTIAL, ELECTRIC FIELD, IN PROBES

      ! Ran into an error restarting from a checkpoint file with zero probes (9/22/23). If there are no probes in the
      !   original simulation the probe data files are length zero. This causes an error when trying to read the first
      !   1:N_of_saved_records data files.
      ! The following line of code checks for a component of this error and skips the reading of the data files if
      !   there are no probes in the continued run:
      !
      !  **** IF ((exists).AND.(N_of_probes.GT.0)) THEN ***
      !
      ! This is not a fool proof fix. If you initially ran with no probes, but now run with probes you will still try
      !   read the first 1:N_of_saved_records lines of the probe data files. This will cause an error.
      ! A way to fix this would be to write the probe data files with a header line that contains the number of probes,
      !   or to write the prior number of probes into the snapshot file, thus allowing us to check how many probes
      !   were used in the prior simulation. If the number of probes was noon-zero, we could proceed with reading the
      !   probe data files. If the number of probes was zero, we could skip the reading of the probe data files and append
      !   any probe data from the continued simulation onto blank files.

! electrons
      INQUIRE (FILE = 'dim_Ne_vst.dat', EXIST = exists)
      IF ((exists).AND.(N_of_probes.GT.0)) THEN
         OPEN (50, FILE = 'dim_Ne_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,999(1x,e12.5))') r_dummy !, arr_dummy(1:N_of_probes)
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! ions
      INQUIRE (FILE = 'dim_Ni_vst.dat', EXIST = exists)
      IF ((exists).AND.(N_of_probes.GT.0)) THEN
         OPEN (50, FILE = 'dim_Ni_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,999(1x,e12.5))') r_dummy !, arr_dummy(1:N_of_probes)
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! potential
      INQUIRE (FILE = 'dim_F_vst.dat', EXIST = exists)
      IF ((exists).AND.(N_of_probes.GT.0)) THEN
         OPEN (50, FILE = 'dim_F_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,999(1x,e12.5))') r_dummy !, arr_dummy(1:N_of_probes)
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! electric field
      INQUIRE (FILE = 'dim_Ex_vst.dat', EXIST = exists)
      IF ((exists).AND.(N_of_probes.GT.0)) THEN
         OPEN (50, FILE = 'dim_Ex_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,999(1x,e12.5))') r_dummy !, arr_dummy(1:N_of_probes)
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! COLLISION FREQUENCIES

! frequencies of electron-neutral collisions
      INQUIRE (FILE = 'dim_fen_collisions_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_fen_collisions_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,5(2x,e12.5))') r_dummy !, r_dummy, r_dummy, r_dummy, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! frequencies of ion-neutral collisions
      INQUIRE (FILE = 'dim_fin_collisions_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_fin_collisions_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,3(2x,e12.5))') r_dummy !, r_dummy, r_dummy, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! AT THE WALLS

! number of particles that hit the left wall: electrons / ions
      INQUIRE (FILE = 'dim_n_leftwall_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_n_leftwall_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,8(2x,e12.5))') r_dummy !, r_dummy, r_dummy, r_dummy, r_dummy, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! number of particles that hit the right wall: electrons / ions
      INQUIRE (FILE = 'dim_n_rightwall_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_n_rightwall_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,8(2x,e12.5))') r_dummy !, r_dummy, r_dummy, r_dummy, r_dummy, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! average energy of particles that hit the left wall: electrons / ions
      INQUIRE (FILE = 'dim_avgenergy_leftwall_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_avgenergy_leftwall_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,2(2x,e12.5))') r_dummy !, r_dummy, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! average energy of particles that hit the right wall: electrons / ions
      INQUIRE (FILE = 'dim_avgenergy_rightwall_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_avgenergy_rightwall_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,2(2x,e12.5))') r_dummy !, r_dummy, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! average energy of particles that are emitted from the left wall: electrons / ions
      INQUIRE (FILE = 'dim_avgenergy_leftemit_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_avgenergy_leftemit_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,3(2x,e12.5))') r_dummy !, r_dummy, r_dummy, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! average energy of particles that are emitted from the right wall: electrons / ions
      INQUIRE (FILE = 'dim_avgenergy_rightemit_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_avgenergy_rightemit_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,3(2x,e12.5))') r_dummy !, r_dummy, r_dummy, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

      ! emission coefficient at the left wall: electrons / ions
      INQUIRE (FILE = 'dim_gamma_left_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_gamma_left_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,e12.5,5(2x,f7.4))') r_dummy !, r_dummy, r_dummy, r_dummy, r_dummy, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! emission coefficient at the right wall: electrons / ions
      INQUIRE (FILE = 'dim_gamma_right_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_gamma_right_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,e12.5,5(2x,f7.4))') r_dummy !, r_dummy, r_dummy, r_dummy, r_dummy, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! number of electrons that hit the left wall, components: secondary / collided / unrecognized
      INQUIRE (FILE = 'dim_n_leftwall_comp_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_n_leftwall_comp_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,3(2x,e12.5))') r_dummy !, r_dummy, r_dummy, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! average energy of el-ns that hit the left wall, components: secondary / collided / unrecognized
      INQUIRE (FILE = 'dim_avgenergy_leftwall_comp_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_avgenergy_leftwall_comp_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,3(2x,e12.5))') r_dummy !, r_dummy, r_dummy, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! emission coefficient at the left wall, components: secondary / collided / unrecognized
      INQUIRE (FILE = 'dim_gamma_left_comp_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_gamma_left_comp_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,3(2x,e12.5))') r_dummy !, r_dummy, r_dummy, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! number of electrons that hit the right wall, components: secondary / collided / unrecognized
      INQUIRE (FILE = 'dim_n_rightwall_comp_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_n_rightwall_comp_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,3(2x,e12.5))') r_dummy !, r_dummy, r_dummy, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! average energy of el-ns that hit the right wall, components: secondary / collided / unrecognized
      INQUIRE (FILE = 'dim_avgenergy_rightwall_comp_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_avgenergy_rightwall_comp_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,3(2x,e12.5))') r_dummy !, r_dummy, r_dummy, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

! emission coefficient at the right wall, components: secondary / collided / unrecognized
      INQUIRE (FILE = 'dim_gamma_right_comp_vst.dat', EXIST = exists)
      IF (exists) THEN
         OPEN (50, FILE = 'dim_gamma_right_comp_vst.dat', STATUS = 'OLD')
         DO i = 1, N_of_saved_records
            READ (50, '(2x,f12.5,3(2x,e12.5))') r_dummy !, r_dummy, r_dummy, r_dummy
         END DO
         ENDFILE 50
         CLOSE (50, STATUS = 'KEEP')
      END IF

   END IF

END SUBROUTINE PREPARE_TIME_DEPENDENCE_DATAFILES

!---------------------------------
!
SUBROUTINE DO_DIAGNOSTICS_STEP_1

   USE ParallelOperationValues
   USE CurrentProblemValues
   USE Diagnostics

   USE MCCollisions, ONLY : e_n_1_count, e_n_2_count, e_n_3_count, e_t_4_count, e_n_5_count, i_n_1_count, i_n_2_count, i_t_3_count

   USE SEEmission, ONLY :   electron_reflux_count,       electron_reflux_energy, &
   & electron_reemit_count,       electron_reemit_energy, &
   & ion_reflux_count,            ion_reflux_energy, &
   & ion_reemit_count,            ion_reemit_energy, &
   & ion_left_reflect_count,      ion_left_reflect_energy, &
   & ion_right_reflect_count,     ion_right_reflect_energy, &
   & see_left_elastic_count,      see_left_elastic_energy, &
   & see_left_inelastic_count,    see_left_inelastic_energy, &
   & see_left_true_count,         see_left_true_energy, &
   & see_right_elastic_count,     see_right_elastic_energy, &
   & see_right_inelastic_count,   see_right_inelastic_energy, &
   & see_right_true_count,        see_right_true_energy, &
   & PlasmaSourceFlag,            Ion_interac_model, &
   & prie_left_from_right_count,  prie_left_from_right_energy, &
   & prie_left_after_coll_count,  prie_left_after_coll_energy, &
   & prie_right_from_left_count,  prie_right_from_left_energy, &
   & prie_right_after_coll_count, prie_right_after_coll_energy, &
   & sece_left_from_right_count,  sece_left_after_coll_count, &
   & sece_right_from_left_count,  sece_right_after_coll_count

   ! use Snapshots, only :  irwedf, ilwedf  ! TEMP-EDF-CLEAR

   IMPLICIT NONE

   INTEGER n       ! index for the bufers of the potential in the midplane and at x=0
   INTEGER s       ! species
   INTEGER k, tag, j       ! index of particle

   INTEGER left_node, right_node, node  !
   REAL(8) x                      ! temporary variables, used for convenience
   REAL(8) vx, vy, vz             !
   REAL(8) dn_left, dn_right, dq_left, dq_right                         ! charge density at neighboring nodes
   REAL(8) jx_left, jy_left, jz_left, jx_right, jy_right, jz_right      ! current elements at neighboring nodes

   IF (Rank_of_process.EQ.0) THEN          ! the server process is responsible for updating the buffers of the potential
      DO n = length_of_fbufer,2,-1
         fmid_bufer( n) = fmid_bufer( n-1)
         fwall_bufer(n) = fwall_bufer(n-1)
      END DO
      fmid_bufer( 1) = F(N_cells/2)
      fwall_bufer(1) = F(0)
   END IF
   ! if (T_cntr.EQ.WriteStart_step) then   ! TEMP-EDF-CLEAR
   !    irwedf = 0.0_8                     ! TEMP-EDF-CLEAR
   !    ilwedf = 0.0_8                     ! TEMP-EDF-CLEAR
   ! end if                                ! TEMP-EDF-CLEAR

! exit if the current time layer is out of the diagnostics time window
   IF (T_cntr.LT.Start_diag_Tcntr) RETURN

   IF (Rank_of_process.GT.0) THEN  ! accumulate data if the process is a client

! cycle over all species
      DO s = 1, N_spec
! cycle over all particles of block, if the block must be accounted
         DO k = 1, N_part(s)
! initiate the temporary variables
            x          =  X_of_spec(s)%part(k)
            left_node  = INT(x)
            IF (left_node.EQ.N_cells) THEN         ! if particle is at the right node exactly (but not collided yet)
               left_node = left_node - 1
            END IF
            right_node = left_node + 1
            vx         = VX_of_spec(s)%part(k)
            vy         = VY_of_spec(s)%part(k)
            vz         = VZ_of_spec(s)%part(k)
            tag        = Tag_of_spec(s)%part(k)

! collect the particle data in the diagnostics arrays
            VX2_cell(right_node, s) = VX2_cell(right_node, s) + vx * vx
            VY2_cell(right_node, s) = VY2_cell(right_node, s) + vy * vy
            VZ2_cell(right_node, s) = VZ2_cell(right_node, s) + vz * vz

            VX_cell(right_node, s) = VX_cell(right_node, s) + vx
            VY_cell(right_node, s) = VY_cell(right_node, s) + vy
            VZ_cell(right_node, s) = VZ_cell(right_node, s) + vz

            Npart_cell(right_node, s) = Npart_cell(right_node, s) + 1

            ! Charge density at neighboring nodes
            dn_right = x - left_node
            dn_left  = right_node - x
            dq_left  = Qs(s) * dn_left
            dq_right = Qs(s) * dn_right

            ! Current elements at neighboring nodes
            jx_left = dq_left * vx
            jy_left = dq_left * vy
            jz_left = dq_left * vz
            jx_right = dq_right * vx
            jy_right = dq_right * vy
            jz_right = dq_right * vz

            QVX_mesh(left_node) = QVX_mesh(left_node) + jx_left
            QVY_mesh(left_node) = QVY_mesh(left_node) + jy_left
            QVZ_mesh(left_node) = QVZ_mesh(left_node) + jz_left

            QVX_mesh(right_node) = QVX_mesh(right_node) + jx_right
            QVY_mesh(right_node) = QVY_mesh(right_node) + jy_right
            QVZ_mesh(right_node) = QVZ_mesh(right_node) + jz_right

            ! Check units on next line before implementing
            ! Rate_energy_heat(s) = Rate_energy_heat(s) + (jx_left * J_scl_Am2) * (EX(left_node) * E_scl_Vm) + &
            !                       & (jx_right * J_scl_Am2) * (EX(right_node) * * E_scl_Vm)

! **** fluxes are updated in the pusher and various emission functions **********
!           NVX_mesh(left_node,s) = NVX_mesh(left_node,s) + dn_left * vx
!           NVX_mesh(right_node,s) = NVX_mesh(right_node,s) + dn_right * vx
! calculate x-fluxes so that the number of particles is conserved,
! right now only works in interior nodes and w/o periodic b.c.:
!           if (tag .ge. 1) then
!             do j = 1, tag
!               node = left_node - j + 1
!               if (node.gt.0) NVX_mesh(node,s) = NVX_mesh(node,s) + 1.  ! crossed in positive direction
!             end do
!           else
!           end if
!           if (tag .le. -1) then
!             do j = 1, -tag
!               node = right_node + j - 1
!               if (node.lt.N_cells) NVX_mesh(node,s) = NVX_mesh(node,s) - 1. ! crossed in negative direction
!             end do
!           else
!           end if

            IF (T_cntr.EQ.Finish_diag_Tcntr) THEN

               Energy_pot(s) = Energy_pot(s) + dq_left * F(left_node) + dq_right * F(right_node)

               Avg_kin_energy_x(s) = Avg_kin_energy_x(s) + vx * vx
               Avg_kin_energy_y(s) = Avg_kin_energy_y(s) + vy * vy
               Avg_kin_energy_z(s) = Avg_kin_energy_z(s) + vz * vz

               VY_recent(s) = VY_recent(s) + vy
               VZ_recent(s) = VZ_recent(s) + vz

            END IF

         END DO
      END DO

   END IF

! if diagnostics in the current time window is about to be completed ...
   IF (T_cntr.EQ.Finish_diag_Tcntr) THEN

! perform output of diagnostic values
      CALL PROCESS_DIAGNOSTIC_DATA  !###

! create the snapshot (if necessary)
!     CALL CREATE_SNAPSHOT

   END IF

END SUBROUTINE DO_DIAGNOSTICS_STEP_1

!---------------------------------
!
SUBROUTINE DO_DIAGNOSTICS_STEP_2

   USE ParallelOperationValues
   USE CurrentProblemValues
   USE Diagnostics

   USE MCCollisions, ONLY : e_n_1_count, e_n_2_count, e_n_3_count, e_t_4_count, e_n_5_count, i_n_1_count, i_n_2_count, i_t_3_count

   USE SEEmission, ONLY :   electron_reflux_count,       electron_reflux_energy, &
   & electron_reemit_count,       electron_reemit_energy, &
   & ion_reflux_count,            ion_reflux_energy, &
   & ion_reemit_count,            ion_reemit_energy, &
   & ion_left_reflect_count,      ion_left_reflect_energy, &
   & ion_right_reflect_count,     ion_right_reflect_energy, &
   & see_left_elastic_count,      see_left_elastic_energy, &
   & see_left_inelastic_count,    see_left_inelastic_energy, &
   & see_left_true_count,         see_left_true_energy, &
   & see_right_elastic_count,     see_right_elastic_energy, &
   & see_right_inelastic_count,   see_right_inelastic_energy, &
   & see_right_true_count,        see_right_true_energy, &
   & PlasmaSourceFlag,            Ion_interac_model, &
   & prie_left_from_right_count,  prie_left_from_right_energy, &
   & prie_left_after_coll_count,  prie_left_after_coll_energy, &
   & prie_right_from_left_count,  prie_right_from_left_energy, &
   & prie_right_after_coll_count, prie_right_after_coll_energy, &
   & sece_left_from_right_count,  sece_left_after_coll_count, &
   & sece_right_from_left_count,  sece_right_after_coll_count

   USE IonInducedSEEmission, ONLY : ionsee_left_count,  ionsee_right_count, &
   & ionsee_left_energy, ionsee_right_energy

   IMPLICIT NONE

! if diagnostics in the current time window is completed ...
   IF (T_cntr.EQ.Finish_diag_Tcntr) THEN

! modify the boundaries of the diagnostics time window
      Start_diag_Tcntr = Start_diag_Tcntr + WriteOut_step
      Finish_diag_Tcntr = Start_diag_Tcntr + WriteAvg_step - 1

! clear the diagnostics arrays
      VX2_cell = 0.0_8
      VY2_cell = 0.0_8
      VZ2_cell = 0.0_8

      VX_cell = 0.0_8
      VY_cell = 0.0_8
      VZ_cell = 0.0_8

      Npart_cell = 0

      QF_mesh = 0.0_8

      QVX_mesh = 0.0_8
      QVY_mesh = 0.0_8
      QVZ_mesh = 0.0_8

      N_new_cell = 0
      P_heat_cell = 0.0_8
      NVX_mesh = 0.0_8

! clear the (energy) quantities accumulated only during the last timestep of the averaging period
      Energy_pot       = 0.0_8    !
      Avg_kin_energy_x = 0.0_8    !
      Avg_kin_energy_y = 0.0_8    ! these values are accumulated in Do_Diagnostics
      Avg_kin_energy_z = 0.0_8    !
      VY_recent        = 0.0_8    !
      VZ_recent        = 0.0_8    !

! clear the energy rates for collisions and heating
      Rate_energy_coll = 0.0_8           ! array
      Rate_energy_heat = 0.0_8

! clear the numbers/energies accumulated on different walls in (final/pre push)
      Rate_number_leftwall  = 0          ! array
      Rate_number_rightwall = 0          ! array
      Rate_number_leftemit  = 0          ! array
      Rate_number_rightemit = 0          ! array
      Rate_energy_leftwall  = 0.0_8      ! array
      Rate_energy_rightwall = 0.0_8      ! array
      Rate_energy_leftemit  = 0.0_8      ! array
      Rate_energy_rightemit = 0.0_8      ! array

      e_n_1_count = 0; e_n_2_count = 0; e_n_3_count = 0; i_n_1_count = 0; i_n_2_count = 0
      e_t_4_count = 0; i_t_3_count = 0
      e_n_5_count = 0

      electron_reflux_count  = 0;  electron_reflux_energy = 0.0_8
      electron_reemit_count  = 0;  electron_reemit_energy = 0.0_8
      ion_reflux_count       = 0;  ion_reflux_energy      = 0.0_8
      ion_reemit_count       = 0;  ion_reemit_energy      = 0.0_8

      see_left_elastic_count   = 0; see_left_elastic_energy   = 0.0_8
      see_left_inelastic_count = 0; see_left_inelastic_energy = 0.0_8
      see_left_true_count      = 0; see_left_true_energy      = 0.0_8

      see_right_elastic_count   = 0; see_right_elastic_energy   = 0.0_8
      see_right_inelastic_count = 0; see_right_inelastic_energy = 0.0_8
      see_right_true_count      = 0; see_right_true_energy      = 0.0_8

      prie_left_from_right_count = 0; prie_left_from_right_energy = 0.0_8
      prie_left_after_coll_count = 0; prie_left_after_coll_energy = 0.0_8

      prie_right_from_left_count  = 0; prie_right_from_left_energy  = 0.0_8
      prie_right_after_coll_count = 0; prie_right_after_coll_energy = 0.0_8

      sece_left_from_right_count = 0; sece_left_after_coll_count = 0
      sece_right_from_left_count = 0; sece_right_after_coll_count = 0

      ionsee_left_count  = 0; ionsee_left_energy = 0.0_8
      ionsee_right_count = 0; ionsee_right_energy = 0.0_8

!     CALL SAVE_CHECKPOINT

   END IF

END SUBROUTINE DO_DIAGNOSTICS_STEP_2

!------------------------------
SUBROUTINE PROCESS_DIAGNOSTIC_DATA
!! We enter this subroutine if we are at a Finish_diag_Tcntr time step. Nodes send all
!! diagnostic data to the server and the server processes it and saves it to files.

   use mpi
   USE ParallelOperationValues
   USE Diagnostics
   USE CurrentProblemValues
   USE MCCollisions
   USE SEEmission
   USE IonInducedSEEmission

   IMPLICIT NONE

!  INCLUDE 'mpif.h'

   INTEGER i     ! node
   INTEGER s     ! species
!  INTEGER n     ! process

   REAL(8) old_Energy_full_eV
   REAL(8) old_Energy_kin_eV(1 : N_spec) ! was (1:2)
   REAL(8) old_Energy_pot_eV(1:N_spec)
   REAL(8) old_Energy_heat_eV(1:N_spec)
   REAL(8) t_s, old_Sigma_electrode, dSigma, j_Am2, P_Wm2 !** for capacitive discharge (BC_flag = 4)

   REAL(8) time_ns
   REAL(8) GetAverageValue
   REAL(8) GetAverageValueInt
   REAL(8) GetAverageFromBufer
   REAL(8) U_ext_vst, delta_t_plot, U_midstep_V

   INTEGER prie_not_recognized_count        ! number of primary electrons, which were not recognized as former secondary (tag = -+1) or collided (tag=2)
   REAL(8) prie_not_recognized_energy       ! energy -"-
   INTEGER sece_not_recognized_count        ! number of secondary electrons, produced by the not recognized primary electrons

   INTEGER ierr
   INTEGER ibufer(1:42)
   REAL(8) dbufer(1:42)
   INTEGER ibufer2(1:42)
   REAL(8) dbufer2(1:42)
!  INTEGER stattus(MPI_STATUS_SIZE)
!  INTEGER tag, source, dest
!  LOGICAL flag

   REAL(8) temp_array(1:999)     ! to gather values in probes before saving to files

   REAL(8) Factor_Qwall_Cm2

   REAL(8) J_external_A, J_Ohm_A

   delta_t_plot = delta_t_s * dble(WriteOut_step)

   ibufer2 = 0
   dbufer2 = 0.0_8

   IF (Rank_of_process.NE.0) THEN                                                                 ! the client >>>
! prepare integer buffer for transmission
      ibufer( 1) = N_part(1)
      ibufer( 2) = N_part(2)
      ibufer( 3) = Q_left
      ibufer( 4) = Q_right
      ibufer( 5) = see_left_elastic_count
      ibufer( 6) = see_left_inelastic_count
      ibufer( 7) = see_left_true_count
      ibufer( 8) = see_right_elastic_count
      ibufer( 9) = see_right_inelastic_count
      ibufer(10) = see_right_true_count
      ibufer(11) = Rate_number_leftwall(1)
      ibufer(12) = Rate_number_leftwall(2)
      ibufer(13) = Rate_number_leftemit(1)
      ibufer(14) = Rate_number_leftemit(2)
      ibufer(15) = Rate_number_rightwall(1)
      ibufer(16) = Rate_number_rightwall(2)
      ibufer(17) = Rate_number_rightemit(1)
      ibufer(18) = Rate_number_rightemit(2)
      ibufer(19) = electron_reflux_count
      ibufer(20) = electron_reemit_count
      ibufer(21) = ion_reflux_count
      ibufer(22) = ion_reemit_count
      ibufer(23) = ion_left_reflect_count
      ibufer(24) = ion_right_reflect_count
      ibufer(25) = e_n_1_count
      ibufer(26) = e_n_2_count
      ibufer(27) = e_n_3_count
      ibufer(28) = e_t_4_count
      ibufer(29) = e_n_5_count
      ibufer(30) = i_n_1_count
      ibufer(31) = i_n_2_count
      ibufer(32) = i_t_3_count
      ibufer(33) = prie_left_from_right_count
      ibufer(34) = prie_left_after_coll_count
      ibufer(35) = prie_right_from_left_count
      ibufer(36) = prie_right_after_coll_count
      ibufer(37) = sece_left_from_right_count
      ibufer(38) = sece_left_after_coll_count
      ibufer(39) = sece_right_from_left_count
      ibufer(40) = sece_right_after_coll_count
      ibufer(41) = ionsee_left_count
      ibufer(42) = ionsee_right_count

! prepare the real(8) buffer for transmission

      dbufer( 1) = VY_recent(1)
      dbufer( 2) = VY_recent(2)
      dbufer( 3) = VZ_recent(1)
      dbufer( 4) = VZ_recent(2)
      dbufer( 5) = see_left_elastic_energy
      dbufer( 6) = see_left_inelastic_energy
      dbufer( 7) = see_left_true_energy
      dbufer( 8) = see_right_elastic_energy
      dbufer( 9) = see_right_inelastic_energy
      dbufer(10) = see_right_true_energy
      dbufer(11) = Rate_energy_leftwall(1)
      dbufer(12) = Rate_energy_leftwall(2)
      dbufer(13) = Rate_energy_leftemit(1)
      dbufer(14) = Rate_energy_leftemit(2)
      dbufer(15) = Rate_energy_rightwall(1)
      dbufer(16) = Rate_energy_rightwall(2)
      dbufer(17) = Rate_energy_rightemit(1)
      dbufer(18) = Rate_energy_rightemit(2)
      dbufer(19) = electron_reflux_energy
      dbufer(20) = electron_reemit_energy
      dbufer(21) = ion_reflux_energy
      dbufer(22) = ion_reemit_energy
      dbufer(23) = ion_left_reflect_energy
      dbufer(24) = ion_right_reflect_energy
      dbufer(25) = Rate_energy_heat(1)
      dbufer(26) = Rate_energy_heat(2)
      dbufer(27) = Rate_energy_coll(1)
      dbufer(28) = Rate_energy_coll(2)
      dbufer(29) = Energy_pot(1)
      dbufer(30) = Energy_pot(2)
      dbufer(31) = Avg_kin_energy_x(1)
      dbufer(32) = Avg_kin_energy_x(2)
      dbufer(33) = Avg_kin_energy_y(1)
      dbufer(34) = Avg_kin_energy_y(2)
      dbufer(35) = Avg_kin_energy_z(1)
      dbufer(36) = Avg_kin_energy_z(2)
      dbufer(37) = prie_left_from_right_energy
      dbufer(38) = prie_left_after_coll_energy
      dbufer(39) = prie_right_from_left_energy
      dbufer(40) = prie_right_after_coll_energy
      dbufer(41) = ionsee_left_energy
      dbufer(42) = ionsee_right_energy

! transmit
      CALL MPI_REDUCE(ibufer, ibufer2, 42, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

      CALL MPI_REDUCE(dbufer, dbufer2, 41, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

! quit
      RETURN

   ELSE                                                                                           ! the server >>>

! clear the buffers
      ibufer = 0
      dbufer = 0.0_8

! receive
!print '(2x,"Process ",i3," will receive IBUFER from all processes")', Rank_of_process
      CALL MPI_REDUCE(ibufer2, ibufer, 42, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

!print '(2x,"Process ",i3," will receive DBUFER from all processes")', Rank_of_process
      CALL MPI_REDUCE(dbufer2, dbufer, 42, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

! restore integer data from the buffer

      N_part(1)                   = ibufer( 1)
      N_part(2)                   = ibufer( 2)
      Q_left                      = ibufer( 3)
      Q_right                     = ibufer( 4)
      see_left_elastic_count      = ibufer( 5)
      see_left_inelastic_count    = ibufer( 6)
      see_left_true_count         = ibufer( 7)
      see_right_elastic_count     = ibufer( 8)
      see_right_inelastic_count   = ibufer( 9)
      see_right_true_count        = ibufer(10)
      Rate_number_leftwall(1)     = ibufer(11)
      Rate_number_leftwall(2)     = ibufer(12)
      Rate_number_leftemit(1)     = ibufer(13)
      Rate_number_leftemit(2)     = ibufer(14)
      Rate_number_rightwall(1)    = ibufer(15)
      Rate_number_rightwall(2)    = ibufer(16)
      Rate_number_rightemit(1)    = ibufer(17)
      Rate_number_rightemit(2)    = ibufer(18)
      electron_reflux_count       = ibufer(19)
      electron_reemit_count       = ibufer(20)
      ion_reflux_count            = ibufer(21)
      ion_reemit_count            = ibufer(22)
      ion_left_reflect_count      = ibufer(23)
      ion_right_reflect_count     = ibufer(24)
      e_n_1_count                 = ibufer(25)
      e_n_2_count                 = ibufer(26)
      e_n_3_count                 = ibufer(27)
      e_t_4_count                 = ibufer(28)
      e_n_5_count                 = ibufer(29)
      i_n_1_count                 = ibufer(30)
      i_n_2_count                 = ibufer(31)
      i_t_3_count                 = ibufer(32)
      prie_left_from_right_count  = ibufer(33)
      prie_left_after_coll_count  = ibufer(34)
      prie_right_from_left_count  = ibufer(35)
      prie_right_after_coll_count = ibufer(36)
      sece_left_from_right_count  = ibufer(37)
      sece_left_after_coll_count  = ibufer(38)
      sece_right_from_left_count  = ibufer(39)
      sece_right_after_coll_count = ibufer(40)
      ionsee_left_count           = ibufer(41)
      ionsee_right_count          = ibufer(42)

! restore real(8) data from the buffer

      VY_recent(1)                 = dbufer( 1)
      VY_recent(2)                 = dbufer( 2)
      VZ_recent(1)                 = dbufer( 3)
      VZ_recent(2)                 = dbufer( 4)
      see_left_elastic_energy      = dbufer( 5)
      see_left_inelastic_energy    = dbufer( 6)
      see_left_true_energy         = dbufer( 7)
      see_right_elastic_energy     = dbufer( 8)
      see_right_inelastic_energy   = dbufer( 9)
      see_right_true_energy        = dbufer(10)
      Rate_energy_leftwall(1)      = dbufer(11)
      Rate_energy_leftwall(2)      = dbufer(12)
      Rate_energy_leftemit(1)      = dbufer(13)
      Rate_energy_leftemit(2)      = dbufer(14)
      Rate_energy_rightwall(1)     = dbufer(15)
      Rate_energy_rightwall(2)     = dbufer(16)
      Rate_energy_rightemit(1)     = dbufer(17)
      Rate_energy_rightemit(2)     = dbufer(18)
      electron_reflux_energy       = dbufer(19)
      electron_reemit_energy       = dbufer(20)
      ion_reflux_energy            = dbufer(21)
      ion_reemit_energy            = dbufer(22)
      ion_left_reflect_energy      = dbufer(23)
      ion_right_reflect_energy     = dbufer(24)
      Rate_energy_heat(1)          = dbufer(25)
      Rate_energy_heat(2)          = dbufer(26)
      Rate_energy_coll(1)          = dbufer(27)
      Rate_energy_coll(2)          = dbufer(28)
      Energy_pot(1)                = dbufer(29)   !    sum of dQ * F, electrons
      Energy_pot(2)                = dbufer(30)   !    sum of dQ * F, ions
      Avg_kin_energy_x(1)          = dbufer(31)   !    sum of VX**2, electrons
      Avg_kin_energy_x(2)          = dbufer(32)   !    sum of VX**2, ions
      Avg_kin_energy_y(1)          = dbufer(33)   !    sum of VY**2, electrons
      Avg_kin_energy_y(2)          = dbufer(34)   !    sum of VY**2, ions
      Avg_kin_energy_z(1)          = dbufer(35)   !    sum of VZ**2, electrons
      Avg_kin_energy_z(2)          = dbufer(36)   !    sum of VZ**2, ions
      prie_left_from_right_energy  = dbufer(37)
      prie_left_after_coll_energy  = dbufer(38)
      prie_right_from_left_energy  = dbufer(39)
      prie_right_after_coll_energy = dbufer(40)
      ionsee_left_energy           = dbufer(41)
      ionsee_right_energy          = dbufer(42)

   END IF

! remember the old values
   old_Energy_full_eV = Energy_full_eV
   old_Energy_kin_eV  = Energy_kin_eV
   old_Energy_pot_eV  = Energy_pot_eV
   old_Energy_heat_eV = Energy_heat_eV
   old_Sigma_electrode = Sigma_electrode

   t_s = dble(T_cntr) * delta_t_s
   U_app_V =  U_ext_vst(t_s) * U_scl_V
   U_midstep_V = U_ext_vst(t_s - 0.5 * delta_t_plot) * U_scl_V
   if (BC_flag.eq.4) then
      U_cap_V =  U_app_V - ( F(0) * U_scl_V )
      Sigma_electrode = eps_0_Fm * eps_layer * (U_cap_V /d_m)
      dSigma = Sigma_electrode - old_Sigma_electrode
      Enr_consumed_Jm2 = Enr_consumed_Jm2 + (U_midstep_V * dSigma)
      j_Am2 = dSigma / (delta_t_plot)
      P_Wm2 = j_Am2 * U_midstep_V
      if (old_Sigma_electrode .eq. 0.0_8) then
         Enr_consumed_Jm2 = 0.0_8 !to start the counter somehow
         j_Am2 = 0.
         P_Wm2 = 0.
      else
      endif
   else
   end if

! calculate kinetic energy for electrons and ions: total / average per particle (APP) / APP, X-motion / APP, Y-motion / APP, Z-motion
   DO s = 1, N_spec
      Avg_kin_energy_x_eV(s) = Factor_energy_eV * MS(s) * Avg_kin_energy_x(s)
      Avg_kin_energy_y_eV(s) = Factor_energy_eV * MS(s) * Avg_kin_energy_y(s)
      Avg_kin_energy_z_eV(s) = Factor_energy_eV * MS(s) * Avg_kin_energy_z(s)

      Energy_kin_eV(s)       = (Avg_kin_energy_x_eV(s) + Avg_kin_energy_y_eV(s) + Avg_kin_energy_z_eV(s)) * N_in_macro

      Avg_kin_energy_x_eV(s) = GetAverageValue(Avg_kin_energy_x_eV(s), N_part(s))
      Avg_kin_energy_y_eV(s) = GetAverageValue(Avg_kin_energy_y_eV(s), N_part(s))
      Avg_kin_energy_z_eV(s) = GetAverageValue(Avg_kin_energy_z_eV(s), N_part(s))

      Avg_kin_energy_eV(s)   = Avg_kin_energy_x_eV(s) + Avg_kin_energy_y_eV(s) + Avg_kin_energy_z_eV(s)
   END DO

! calculate the potential energy
   DO s = 1, N_spec
      Energy_pot_eV(s) = Factor_energy_pot_eV * (Energy_pot(s) + Q_left * F(0) + Q_right * F(N_cells))
   END DO

! a correction for the walls with given potential; does not include energy stored in dielectric layer
   IF (BC_flag.EQ.0 .or. BC_flag.eq.4) THEN
      DO s = 1, N_spec
         Energy_pot_eV(s) = Energy_pot_eV(s) + 0.5_8 * eps_0_Fm * F(0) * EX(0) * E_scl_Vm / e_Cl
      END DO
   END IF

! new full energy of system
   Energy_full_eV = Energy_kin_eV(1) + Energy_kin_eV(2) + Energy_pot_eV(1) + Energy_pot_eV(2)

! rate of change of full energy, kinetic energy, potential energy
   ! If we are at the first Finish_diag_Tcntr step
   IF (T_cntr.EQ.(WriteStart_step + WriteAvg_step - 1)) THEN   ! For the first time ALL RATES should be set to zero,
      Init_energy_full_eV    = Energy_full_eV                  ! because we initiate the system's full energy here -
      Rate_energy_full_eVns1 = 0.0_8                           ! i.e. create the start point for the energy balance -
      Rate_energy_pot_eVns1  = 0.0_8                           ! all losses / gains of energy will be accounted starting this moment...
      Rate_energy_kin_eVns1  = 0.0_8         ! array           ! Otherwise there is always an error in the energy conservation...
      Rate_energy_wall       = 0.0_8   !
      Rate_energy_emit       = 0.0_8   ! arrays
      Rate_energy_coll       = 0.0_8   !
      Energy_heat_eV         = 0.0_8              ! total energy change due to Joule heating
      Rate_energy_heat_eVns1 = 0.0_8              ! rate of change of energy due to Joule heating, dimensional

   ! If we are at any other Finish_diag_Tcntr step
   ELSE
      Rate_energy_full_eVns1 = (Energy_full_eV - old_Energy_full_eV) * Factor_rate_ns1
      Rate_energy_pot_eVns1  = (Energy_pot_eV  - old_Energy_pot_eV)  * Factor_rate_ns1
      Rate_energy_kin_eVns1  = (Energy_kin_eV  - old_Energy_kin_eV)  * Factor_rate_ns1  ! array
      Rate_energy_wall       = -Rate_energy_rightwall - Rate_energy_leftwall   ! arrays
      Rate_energy_emit       =  Rate_energy_rightemit + Rate_energy_leftemit   !
      Energy_heat_eV         =  Energy_heat_eV + Factor_energy_macro_eV * Rate_energy_heat  ! total energy change due to Joule heating
      Rate_energy_heat_eVns1 = (Energy_heat_eV - old_Energy_heat_eV) * Factor_rate_ns1      ! rate of change of energy due to Joule heating, dimensional

   END IF

! separately for electrons and ions:
   DO s = 1, N_spec

      Energy_wall_eV(s) = Energy_wall_eV(s) + Factor_energy_macro_eV * MS(s) * Rate_energy_wall(s)  ! to the walls                        !>>
      Energy_emit_eV(s) = Energy_emit_eV(s) + Factor_energy_macro_eV * MS(s) * Rate_energy_emit(s)  ! from the walls                      !>>
      Energy_coll_eV(s) = Energy_coll_eV(s) + Factor_energy_macro_eV * MS(s) * Rate_energy_coll(s)  ! collisions                          !>>

      Rate_energy_wall_eVns1(s) = Factor_rate_eVns1 * MS(s) * Rate_energy_wall(s)             ! corresponding rates
      Rate_energy_emit_eVns1(s) = Factor_rate_eVns1 * MS(s) * Rate_energy_emit(s)             !
      Rate_energy_coll_eVns1(s) = Factor_rate_eVns1 * MS(s) * Rate_energy_coll(s)             !

      IF (Rate_number_leftwall(s).GT.0) THEN
         Rate_number_leftwall_ns1(s) = Factor_rate_macro_ns1 * Rate_number_leftwall(s)
         Avg_energy_leftwall_eV(s)   = Factor_energy_eV * MS(s) * Rate_energy_leftwall(s) / Rate_number_leftwall(s)                              !<<
         Yield_left(s)               = DBLE(Rate_number_leftemit(s)) / DBLE(Rate_number_leftwall(s))
         IF (s.EQ.1) THEN
            eYield_left_elast   = DBLE(see_left_elastic_count)   / DBLE(Rate_number_leftwall(s))
            eYield_left_inelast = DBLE(see_left_inelastic_count) / DBLE(Rate_number_leftwall(s))
            eYield_left_true    = DBLE(see_left_true_count)      / DBLE(Rate_number_leftwall(s))
         END IF
      ELSE
         Rate_number_leftwall_ns1(s) = 0.0_8
         Avg_energy_leftwall_eV(s)   = 0.0_8
         Yield_left(s)               = 0.0_8
         IF (s.EQ.1) THEN
            eYield_left_elast   = 0.0_8
            eYield_left_inelast = 0.0_8
            eYield_left_true    = 0.0_8
         END IF
      END IF

      IF (Rate_number_rightwall(s).GT.0) THEN
         Rate_number_rightwall_ns1(s) = Factor_rate_macro_ns1 * Rate_number_rightwall(s)
         Avg_energy_rightwall_eV(s)   = Factor_energy_eV * MS(s) * Rate_energy_rightwall(s) / Rate_number_rightwall(s)                           !<<
         Yield_right(s)               = DBLE(Rate_number_rightemit(s)) / DBLE(Rate_number_rightwall(s))
         IF (s.EQ.1) THEN
            eYield_right_elast   = DBLE(see_right_elastic_count)   / DBLE(Rate_number_rightwall(s))
            eYield_right_inelast = DBLE(see_right_inelastic_count) / DBLE(Rate_number_rightwall(s))
            eYield_right_true    = DBLE(see_right_true_count)      / DBLE(Rate_number_rightwall(s))
         END IF
      ELSE
         Rate_number_rightwall_ns1(s) = 0.0_8
         Avg_energy_rightwall_eV(s)   = 0.0_8
         Yield_right(s)               = 0.0_8
         IF (s.EQ.1) THEN
            eYield_right_elast   = 0.0_8
            eYield_right_inelast = 0.0_8
            eYield_right_true    = 0.0_8
         END IF
      END IF

      Rate_number_leftemit_ns1(s) = Factor_rate_macro_ns1 * Rate_number_leftemit(s)
      Avg_energy_leftemit_eV(s)   = Factor_energy_eV * MS(s) * &
      & GetAverageValue(Rate_energy_leftemit(s), Rate_number_leftemit(s))                 !<<

      Rate_number_rightemit_ns1(s) = Factor_rate_macro_ns1 * Rate_number_rightemit(s)
      Avg_energy_rightemit_eV(s)   = Factor_energy_eV * MS(s) * &
      & GetAverageValue(Rate_energy_rightemit(s), Rate_number_rightemit(s))              !<<

   END DO

! calculate the average Y- and Z-velocities for electrons
   VY_e_avg_ms = VY_recent(1) * v_Te_ms / (DBLE(N_box_vel) * DBLE(N_part(1)))
   VZ_e_avg_ms = VZ_recent(1) * v_Te_ms / (DBLE(N_box_vel) * DBLE(N_part(1)))

! calculate the average Y- and Z-electric current density
   JY_avg_Am2 = e_Cl * N_plasma_m3 * (Qs(1) * VY_recent(1) + Qs(2) * VY_recent(2)) * v_Te_ms / &
   & (DBLE(N_box_vel) * DBLE(N_nodes) * DBLE(N_of_particles_cell))
   JZ_avg_Am2 = e_Cl * N_plasma_m3 * (Qs(1) * VZ_recent(1) + Qs(2) * VZ_recent(2)) * v_Te_ms / &
   & (DBLE(N_box_vel) * DBLE(N_nodes) * DBLE(N_of_particles_cell))

! calculate the average electron and ion density (for simplicity we use the streaming values)
   N_avg_m3 = N_part * N_plasma_m3 / (DBLE(N_of_particles_cell) * DBLE(N_nodes))

! factors for frequencies of collisions with neutrals, must be calculated BEFORE the "DoTextOutput" is called !!!
   DO s = 1, N_spec                                                                   ! LOOK AT THE INTEGERS: NOTE, THE SOURCE OF SERIOUS ERROR WAS FOUND HERE:
      IF (N_part(s).GT.0) THEN                                                        ! THE PRODUCT OF TWO INTEGERS "N_part(s) * WriteOut_step" CAN EXCEED
         f_factor(s) = 1.0_8 / (DBLE(N_part(s)) * DBLE(WriteOut_step) * delta_t_s)    ! THE UPPER LIMIT FOR THE 32-BIT INTEGER 2^31-1 AND THE RESULT,
      ELSE                                                                            ! f_factor, BECOMES NEGATIVE AND IRRELEVANT !!!!
         f_factor(s) = 0.0_8                                                          ! NEXT TIME USE 64-BIT INTEGERS OR DIRECTLY CONVER THEM TO FLOATS
      END IF
   END DO

! modify counter to skip extra diagnostics periods between two text outputs
   text_output_counter = text_output_counter + 1
! perform text output
   IF (text_output_counter.GT.TextOut_avg_prds_to_skip) THEN
      CALL DoTextOutput
      text_output_counter = 0
   END IF

   time_ns = dble(T_cntr) * delta_t_s * 1.0e9
!--

   OPEN  (50, FILE = 'dim_fullenergy_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,2x,e12.5)')    time_ns, Energy_full_eV
   CLOSE (50, STATUS = 'KEEP')

   OPEN  (50, FILE = 'dim_potenergy_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,2x,e12.5)')    time_ns, Energy_pot_eV(1) + Energy_pot_eV(2), Energy_pot_eV(1), Energy_pot_eV(2)
   CLOSE (50, STATUS = 'KEEP')

   OPEN  (50, FILE = 'dim_kinenergy_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,3(2x,e12.5))') time_ns, Energy_kin_eV(1) + Energy_kin_eV(2), Energy_kin_eV(1), Energy_kin_eV(2)
   CLOSE (50, STATUS = 'KEEP')

   OPEN  (50, FILE = 'dim_energy_1e_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,4(2x,e12.5))') time_ns, Avg_kin_energy_eV(1), Avg_kin_energy_x_eV(1), &
   & Avg_kin_energy_y_eV(1), Avg_kin_energy_z_eV(1)
   CLOSE (50, STATUS = 'KEEP')

   OPEN  (50, FILE = 'dim_energy_1i_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,4(2x,e12.5))') time_ns, Avg_kin_energy_eV(2), Avg_kin_energy_x_eV(2), &
   & Avg_kin_energy_y_eV(2), Avg_kin_energy_z_eV(2)
   CLOSE (50, STATUS = 'KEEP')

! kinetic energy of particles that have struck the walls:
   OPEN  (50, FILE = 'dim_energy_wall_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,3(2x,e12.5))') time_ns, Energy_wall_eV(1) + Energy_wall_eV(2), Energy_wall_eV(1), Energy_wall_eV(2)
   CLOSE (50, STATUS = 'KEEP')

! kinetic energy of particles that were emitted/reflected from the walls:
   OPEN  (50, FILE = 'dim_energy_emit_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,3(2x,e12.5))') time_ns, Energy_emit_eV(1) + Energy_emit_eV(2), Energy_emit_eV(1), Energy_emit_eV(2)
   CLOSE (50, STATUS = 'KEEP')

! kinetic energy, lost due to collisions with neutrals:
   OPEN  (50, FILE = 'dim_energy_coll_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,3(2x,e12.5))') time_ns, Energy_coll_eV(1) + Energy_coll_eV(2), Energy_coll_eV(1), Energy_coll_eV(2)
   CLOSE (50, STATUS = 'KEEP')

!! energy of Joule heating of plasma by accelerating electric field
!  OPEN  (50, FILE = 'dim_energy_heat_vst.dat', POSITION = 'APPEND')
!  WRITE (50, '(7(2x, e12.5))')  time_ns, Energy_heat_eV, Enr_consumed_Jm2, j_Am2, U_app_V, U_cap_V, P_Wm2
!  CLOSE (50, STATUS = 'KEEP')

! violation of conservation of energy (%)
   OPEN  (50, FILE = 'dim_energy_cons_violat_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,2x,e12.5)')    time_ns, 100.0_8 * (Energy_full_eV - (Energy_wall_eV(1) + Energy_wall_eV(2) + &
   & Energy_emit_eV(1) + Energy_emit_eV(2) + &
   & Energy_coll_eV(1) + Energy_coll_eV(2) + &
   & Energy_heat_eV(1) + Energy_heat_eV(2)) - Init_energy_full_eV) / Energy_full_eV
   CLOSE (50, STATUS = 'KEEP')

! full energy of system (kinetic + potential)
   OPEN  (50, FILE = 'dim_rate_fullenergy_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,2x,e12.5)')    time_ns, Rate_energy_full_eVns1
   CLOSE (50, STATUS = 'KEEP')

! potential energy of system
   OPEN  (50, FILE = 'dim_rate_potenergy_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,2x,e12.5)')    time_ns, Rate_energy_pot_eVns1(1) + Rate_energy_pot_eVns1(2), &
   & Rate_energy_pot_eVns1(1), Rate_energy_pot_eVns1(2)
   CLOSE (50, STATUS = 'KEEP')

! kinetic energy of system: all species / electrons / ions
   OPEN  (50, FILE = 'dim_rate_kinenergy_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,3(2x,e12.5))') time_ns, Rate_energy_kin_eVns1(1)  + Rate_energy_kin_eVns1(2),  &
   & Rate_energy_kin_eVns1(1),  Rate_energy_kin_eVns1(2)
   CLOSE (50, STATUS = 'KEEP')

! kinetic energy of particles that have struck the walls:
   OPEN  (50, FILE = 'dim_rate_energy_wall_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,3(2x,e12.5))') time_ns, Rate_energy_wall_eVns1(1) + Rate_energy_wall_eVns1(2), &
   & Rate_energy_wall_eVns1(1), Rate_energy_wall_eVns1(2)
   CLOSE (50, STATUS = 'KEEP')

! kinetic energy of particles that were emitted/reflected from the walls:
   OPEN  (50, FILE = 'dim_rate_energy_emit_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,3(2x,e12.5))') time_ns, Rate_energy_emit_eVns1(1) + Rate_energy_emit_eVns1(2), &
   & Rate_energy_emit_eVns1(1), Rate_energy_emit_eVns1(2)
   CLOSE (50, STATUS = 'KEEP')

! kinetic energy, lost due to collisions with neutrals:
   OPEN  (50, FILE = 'dim_rate_energy_coll_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,3(2x,e12.5))') time_ns, Rate_energy_coll_eVns1(1) + Rate_energy_coll_eVns1(2), &
   & Rate_energy_coll_eVns1(1), Rate_energy_coll_eVns1(2)
   CLOSE (50, STATUS = 'KEEP')

! energy of Joule heating of plasma by accelerating electric field
   OPEN  (50, FILE = 'dim_rate_energy_heat_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,2x,e12.5)')    time_ns, Rate_energy_heat_eVns1(1) + Rate_energy_heat_eVns1(2), &
   & Rate_energy_heat_eVns1(1), Rate_energy_heat_eVns1(2)
   CLOSE (50, STATUS = 'KEEP')

! velocity of electron Y-flow
   OPEN  (50, FILE = 'dim_vy_e_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,2x,e12.5)')    time_ns, VY_e_avg_ms
   CLOSE (50, STATUS = 'KEEP')

! velocity of electron Z-flow
   OPEN  (50, FILE = 'dim_vz_e_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,2x,e12.5)')    time_ns, VZ_e_avg_ms
   CLOSE (50, STATUS = 'KEEP')

! electric current density in Y-direction
   OPEN  (50, FILE = 'dim_jy_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,2x,e12.5)')    time_ns, JY_avg_Am2
   CLOSE (50, STATUS = 'KEEP')

! electric current density in Z-direction
   OPEN  (50, FILE = 'dim_jz_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,2x,e12.5)')    time_ns, JZ_avg_Am2
   CLOSE (50, STATUS = 'KEEP')

! VALUES IN PROBES >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

   IF (N_of_probes.GT.0) THEN

! electron density
      OPEN  (50, FILE = 'dim_Ne_vst.dat', POSITION = 'APPEND')
      DO i = 1, N_of_probes
         temp_array(i) = Q_strm_spec(probe_node(i), 1) * N_scl_m3
         IF ((probe_node(i).EQ.0).OR.(probe_node(i).EQ.N_cells)) temp_array(i) = 2.0_8 * temp_array(i)
      END DO
      WRITE (50, '(2x,f12.5,999(1x,e12.5))') time_ns, temp_array(1:N_of_probes)
      CLOSE (50, STATUS = 'KEEP')

! ion density
      IF (N_spec.EQ.2) THEN
         OPEN  (50, FILE = 'dim_Ni_vst.dat', POSITION = 'APPEND')
         DO i = 1, N_of_probes
            temp_array(i) = Q_strm_spec(probe_node(i), 2) * N_scl_m3
            IF ((probe_node(i).EQ.0).OR.(probe_node(i).EQ.N_cells)) temp_array(i) = 2.0_8 * temp_array(i)
         END DO
         WRITE (50, '(2x,f12.5,999(1x,e12.5))') time_ns, temp_array(1:N_of_probes)
         CLOSE (50, STATUS = 'KEEP')
      END IF

! potential
      OPEN  (50, FILE = 'dim_F_vst.dat', POSITION = 'APPEND')
      DO i = 1, N_of_probes
         temp_array(i) = F(probe_node(i)) * U_scl_V
      END DO
      WRITE (50, '(2x,f12.5,999(1x,e12.5))') time_ns, temp_array(1:N_of_probes)
      CLOSE (50, STATUS = 'KEEP')

! electric field
      OPEN  (50, FILE = 'dim_Ex_vst.dat', POSITION = 'APPEND')
      DO i = 1, N_of_probes
         temp_array(i) = EX(probe_node(i)) * E_scl_Vm
      END DO
      WRITE (50, '(2x,f12.5,999(1x,e12.5))') time_ns, temp_array(1:N_of_probes)
      CLOSE (50, STATUS = 'KEEP')

   END IF

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

! frequencies of electron-neutral collisions: elastic / excitation-1 / excitation-2 / ionization / turbulence
   OPEN  (50, FILE = 'dim_fen_collisions_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,5(2x,e12.5))') time_ns, e_n_1_count * f_factor(1), e_n_2_count * f_factor(1),  e_n_5_count * f_factor(1), e_n_3_count * f_factor(1), e_t_4_count * f_factor(1)
   CLOSE (50, STATUS = 'KEEP')

! frequencies of ion-neutral collisions: elastic / charge exchange / turbulence
   OPEN  (50, FILE = 'dim_fin_collisions_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,3(2x,e12.5))') time_ns, i_n_1_count * f_factor(2), i_n_2_count * f_factor(2), i_t_3_count * f_factor(2)
   CLOSE (50, STATUS = 'KEEP')

   Factor_Qwall_Cm2 = N_scl_m3 * delta_x_m * e_Cl

! number of particles that hit the left wall: electrons / ions
   OPEN  (50, FILE = 'dim_n_leftwall_vst.dat', POSITION = 'APPEND')

   IF (BC_flag.EQ.1) THEN                ! floating potential
      J_Ohm_A = 0.0_8                    ! external current should be zero
      full_Q_left = DBLE(Q_left)
   ELSE IF (BC_flag.EQ.0) THEN           ! constant potential, R_external = 0
      J_Ohm_A = 0.0_8                    ! cannot use Ohm's law here
   ELSE IF (BC_flag.EQ.2) THEN           ! external circuit
      J_Ohm_A = (U_ext_vst(t_s) - F(0)) * U_scl_V / R_ext_ohm ! Old implementation, which was incorrect for rf sources: (U_ext - F(0)) * U_scl_V / R_ext_ohm
   END IF

   J_external_A = ( full_Q_left - &
   & prev_Q_left + &
   & Qs(1) * (Rate_number_leftemit(1) - Rate_number_leftwall(1)) + &
   & Qs(2) * (Rate_number_leftemit(2) - Rate_number_leftwall(2)) ) * 1.0d-4 * S_electrode_cm2 * e_Cl * 1.0d9 * Factor_rate_macro_ns1

   IF (BC_flag.EQ.0) then
      j_Am2 = J_external_A / S_electrode_cm2 * 1.0d4
      P_Wm2 = j_Am2 * U_midstep_V
      Enr_consumed_Jm2 = Enr_consumed_Jm2 + P_Wm2 * delta_t_plot
      U_cap_V = 0.0_8
   else
   end if

   WRITE (50, '(2x,f12.5,10(2x,e12.5))') &
   & time_ns, &
   & Rate_number_leftwall_ns1(1), &
   & Rate_number_leftwall_ns1(2), &
   & Rate_number_leftemit_ns1(1), &
   & Rate_number_leftemit_ns1(2), &
   & ionsee_left_count * Factor_rate_macro_ns1, &
   & full_Q_left * Factor_Qwall_Cm2, &    ! surface charge after final push via the near-wall volume charge [contains contrib. from ext. circuit]
   & EX(0) * E_scl_Vm, &                  ! electric field obtained after the preliminary push and used during the final push
   & F(0) * U_scl_V, &                    ! electrostatic potential of the wall
   & J_external_A, &                      ! current in the external circuit from conservation of the surface charge
   & J_Ohm_A                              ! Ohm's law

   prev_Q_left = full_Q_left
   CLOSE (50, STATUS = 'KEEP')

! energy of Joule heating of plasma by accelerating electric field
   OPEN  (50, FILE = 'dim_energy_heat_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(7(2x, e12.5))')  time_ns, Energy_heat_eV(1) + Energy_heat_eV(2), Energy_heat_eV(1), Energy_heat_eV(2)!, Enr_consumed_Jm2, j_Am2, U_app_V, U_cap_V, P_Wm2
   CLOSE (50, STATUS = 'KEEP')

! number of particles that hit the right wall: electrons / ions
   OPEN  (50, FILE = 'dim_n_rightwall_vst.dat', POSITION = 'APPEND')

   IF (BC_flag.EQ.1) THEN                ! floating potential
      J_Ohm_A = 0.0_8                    ! external current should be zero
      full_Q_right = DBLE(Q_right)
   ELSE IF (BC_flag.EQ.0) THEN           ! constant potential, R_external = 0
      J_Ohm_A = 0.0_8                    ! cannot use Ohm's law here
   ELSE IF (BC_flag.EQ.2) THEN           ! external circuit
      J_Ohm_A = (U_ext_vst(t_s) - F(0)) * U_scl_V / R_ext_ohm ! Old implementation, which was incorrect for rf sources: (U_ext - F(0)) * U_scl_V / R_ext_ohm
   END IF

   J_external_A = -( full_Q_right - &
   & prev_Q_right + &
   & Qs(1) * (Rate_number_rightemit(1) - Rate_number_rightwall(1)) + &
   & Qs(2) * (Rate_number_rightemit(2) - Rate_number_rightwall(2)) ) * 1.0d-4 * S_electrode_cm2 * e_Cl * 1.0d9 * Factor_rate_macro_ns1

   WRITE (50, '(2x,f12.5,10(2x,e12.5))') &
   & time_ns, &
   & Rate_number_rightwall_ns1(1), &
   & Rate_number_rightwall_ns1(2), &
   & Rate_number_rightemit_ns1(1), &
   & Rate_number_rightemit_ns1(2), &
   & ionsee_right_count * Factor_rate_macro_ns1, &
   & full_Q_right * Factor_Qwall_Cm2, &
   & EX(N_cells) * E_scl_Vm, &
   & F(N_cells) * U_scl_V, &
   & J_external_A, &
   & J_Ohm_A

   prev_Q_right = full_Q_right
   CLOSE (50, STATUS = 'KEEP')

! average energy of particles that hit the left wall: electrons / ions
   OPEN  (50, FILE = 'dim_avgenergy_leftwall_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,2(2x,e12.5))') time_ns, Avg_energy_leftwall_eV(1), Avg_energy_leftwall_eV(2)
   CLOSE (50, STATUS = 'KEEP')

! average energy of particles that hit the right wall: electrons / ions
   OPEN  (50, FILE = 'dim_avgenergy_rightwall_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,2(2x,e12.5))') time_ns, Avg_energy_rightwall_eV(1), Avg_energy_rightwall_eV(2)
   CLOSE (50, STATUS = 'KEEP')

! average energy of particles that are emitted from the left wall: electrons / ions
   OPEN  (50, FILE = 'dim_avgenergy_leftemit_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,3(2x,e12.5))') &
   & time_ns, &
   & Factor_energy_eV * MS(1) * GetAverageValue(Rate_energy_leftemit(1)-ionsee_left_energy, &
   &                                            Rate_number_leftemit(1)-ionsee_left_count), &
   & Avg_energy_leftemit_eV(2), &
   & Factor_energy_eV * MS(1) * GetAverageValue(ionsee_left_energy, ionsee_left_count)
   CLOSE (50, STATUS = 'KEEP')

! average energy of particles that are emitted from the right wall: electrons / ions
   OPEN  (50, FILE = 'dim_avgenergy_rightemit_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,3(2x,e12.5))') &
   & time_ns, &
   & Factor_energy_eV * MS(1) * GetAverageValue(Rate_energy_rightemit(1)-ionsee_right_energy, &
   &                                            Rate_number_rightemit(1)-ionsee_right_count), &
   & Avg_energy_rightemit_eV(2), &
   & Factor_energy_eV * MS(1) * GetAverageValue(ionsee_right_energy, ionsee_right_count)
   CLOSE (50, STATUS = 'KEEP')

! emission coefficient at the left wall: electrons / ions
   OPEN  (50, FILE = 'dim_gamma_left_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,e12.5,5(2x,f7.4))')  time_ns, Yield_left(1),  Yield_left(2),  eYield_left_elast, &
   & eYield_left_inelast,  eYield_left_true
   CLOSE (50, STATUS = 'KEEP')

! emission coefficient at the right wall: electrons / ions
   OPEN  (50, FILE = 'dim_gamma_right_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,e12.5,5(2x,f7.4))')  time_ns, Yield_right(1), Yield_right(2), eYield_right_elast, &
   & eYield_right_inelast, eYield_right_true
   CLOSE (50, STATUS = 'KEEP')

! number of electrons that hit the left wall
   prie_not_recognized_count  = Rate_number_leftwall(1) - prie_left_from_right_count  - prie_left_after_coll_count
   prie_not_recognized_energy = Rate_energy_leftwall(1) - prie_left_from_right_energy - prie_left_after_coll_energy
   sece_not_recognized_count  = Rate_number_leftemit(1) - sece_left_from_right_count  - sece_left_after_coll_count

   OPEN  (50, FILE = 'dim_n_leftwall_comp_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,3(2x,e12.5))') time_ns, prie_left_from_right_count * Factor_rate_macro_ns1, &
   & prie_left_after_coll_count * Factor_rate_macro_ns1, &
   & prie_not_recognized_count  * Factor_rate_macro_ns1
   CLOSE (50, STATUS = 'KEEP')

! average energy of el-ns that hit the left wall
   OPEN  (50, FILE = 'dim_avgenergy_leftwall_comp_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,3(2x,e12.5))') time_ns, &
   & Factor_energy_eV * MS(1) * GetAverageValue(prie_left_from_right_energy, prie_left_from_right_count), &
   & Factor_energy_eV * MS(1) * GetAverageValue(prie_left_after_coll_energy, prie_left_after_coll_count), &
   & Factor_energy_eV * MS(1) * GetAverageValue(prie_not_recognized_energy,  prie_not_recognized_count)
   CLOSE (50, STATUS = 'KEEP')

! emission coefficient at the left wall
   OPEN  (50, FILE = 'dim_gamma_left_comp_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,3(2x,e12.5))') time_ns, GetAverageValueInt(sece_left_from_right_count, prie_left_from_right_count), &
   & GetAverageValueInt(sece_left_after_coll_count, prie_left_after_coll_count), &
   & GetAverageValueInt(sece_not_recognized_count,  prie_not_recognized_count)
   CLOSE (50, STATUS = 'KEEP')

! number of electrons that hit the right wall
   prie_not_recognized_count  = Rate_number_rightwall(1) - prie_right_from_left_count  - prie_right_after_coll_count
   prie_not_recognized_energy = Rate_energy_rightwall(1) - prie_right_from_left_energy - prie_right_after_coll_energy
   sece_not_recognized_count  = Rate_number_rightemit(1) - sece_right_from_left_count  - sece_right_after_coll_count

   OPEN  (50, FILE = 'dim_n_rightwall_comp_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,3(2x,e12.5))') time_ns, prie_right_from_left_count  * Factor_rate_macro_ns1, &
   & prie_right_after_coll_count * Factor_rate_macro_ns1, &
   & prie_not_recognized_count   * Factor_rate_macro_ns1
   CLOSE (50, STATUS = 'KEEP')

! average energy of el-ns that hit the right wall
   OPEN  (50, FILE = 'dim_avgenergy_rightwall_comp_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,3(2x,e12.5))') time_ns, &
   & Factor_energy_eV * MS(1) * GetAverageValue(prie_right_from_left_energy,  prie_right_from_left_count), &
   & Factor_energy_eV * MS(1) * GetAverageValue(prie_right_after_coll_energy, prie_right_after_coll_count), &
   & Factor_energy_eV * MS(1) * GetAverageValue(prie_not_recognized_energy,   prie_not_recognized_count)
   CLOSE (50, STATUS = 'KEEP')

! emission coefficient at the right wall
   OPEN  (50, FILE = 'dim_gamma_right_comp_vst.dat', POSITION = 'APPEND')
   WRITE (50, '(2x,f12.5,3(2x,e12.5))') time_ns, GetAverageValueInt(sece_right_from_left_count,  prie_right_from_left_count), &
   & GetAverageValueInt(sece_right_after_coll_count, prie_right_after_coll_count), &
   & GetAverageValueInt(sece_not_recognized_count,   prie_not_recognized_count)
   CLOSE (50, STATUS = 'KEEP')

! modify the number of saved records
   N_of_saved_records = N_of_saved_records + 1

END SUBROUTINE PROCESS_DIAGNOSTIC_DATA

!------------------------------
!
SUBROUTINE DoTextOutput

   USE Diagnostics

   USE CurrentProblemValues, ONLY :   N_spec, N_part, delta_t_s, T_cntr, MS, T_e_eV !, WriteAvg_step

   USE MCCollisions, ONLY : e_n_1_count, e_n_2_count, e_n_3_count, e_t_4_count, e_n_5_count, i_n_1_count, i_n_2_count, i_t_3_count

   USE SEEmission, ONLY :   electron_reflux_count,     electron_reflux_energy, &
   & electron_reemit_count,     electron_reemit_energy, &
   & ion_reflux_count,          ion_reflux_energy, &
   & ion_reemit_count,          ion_reemit_energy, &
   & ion_left_reflect_count,    ion_left_reflect_energy, &
   & ion_right_reflect_count,   ion_right_reflect_energy, &
   & see_left_elastic_count,    see_left_elastic_energy, &
   & see_left_inelastic_count,  see_left_inelastic_energy, &
   & see_left_true_count,       see_left_true_energy, &
   & see_right_elastic_count,   see_right_elastic_energy, &
   & see_right_inelastic_count, see_right_inelastic_energy, &
   & see_right_true_count,      see_right_true_energy, &
   & PlasmaSourceFlag,          Ion_interac_model

   IMPLICIT NONE

   INTEGER s       ! species
   INTEGER balance_e
   INTEGER balance_i

   REAL(8) GetAverageValue

   PRINT '(/,"Diagnostics report at step : ",i9,", time interval : ",f12.4," ns ... ",f12.4," ns")', &
   & T_cntr, (T_cntr - WriteAvg_step) * delta_t_s * 1.0e9, &
   & T_cntr * delta_t_s * 1.0e9
   PRINT  '("Parts of full system energy [eV]:         pot. = ",e11.4,", pot.e. = ",e11.4,", pot.i. = ",e11.4)', &
   & Energy_pot_eV(1) + Energy_pot_eV(2), &
   & Energy_pot_eV(1), &
   & Energy_pot_eV(2)
   PRINT  '("                                          kin. = ",e11.4,", kin.e. = ",e11.4,", kin.i. = ",e11.4)', &
   & Energy_kin_eV(1) + Energy_kin_eV(2), &
   & Energy_kin_eV(1), &
   & Energy_kin_eV(2)

   PRINT &
   & '("Energy losses/gains [eV]:                 wall = ",e11.4,", emit. = ",e11.4,",   col. = ",e11.4,",  heat. = ",e11.4)', &
   & Energy_wall_eV(1) + Energy_wall_eV(2), &
   & Energy_emit_eV(1) + Energy_emit_eV(2), &
   & Energy_coll_eV(1) + Energy_coll_eV(2), &
   & Energy_heat_eV(1) + Energy_heat_eV(2)

   PRINT  '("Full system energy (kinetic + potential) :       ",e11.4," eV" )', Energy_full_eV

   PRINT  '("Sum of losses/gains + initial energy     :       ",e11.4," eV" )', &
   & Energy_wall_eV(1) + Energy_wall_eV(2) + &
   & Energy_emit_eV(1) + Energy_emit_eV(2) + &
   & Energy_coll_eV(1) + Energy_coll_eV(2) + &
   & Energy_heat_eV(1) + Energy_heat_eV(2) + &
   & Init_energy_full_eV

   PRINT  '("Violation of full energy conservation law:       ",f11.1," %" )', &
   & 100.0_8 * (Energy_full_eV - (Energy_wall_eV(1) + Energy_wall_eV(2) + &
   & Energy_emit_eV(1) + Energy_emit_eV(2) + &
   & Energy_coll_eV(1) + Energy_coll_eV(2) + &
   & Energy_heat_eV(1) + Energy_heat_eV(2)) - Init_energy_full_eV) / Energy_full_eV

   PRINT &
   & '("Energy rates [eV/ns]: full = ",e11.4,", wall = ",e11.4,", emit. = ",e11.4,",  coll. = ",e11.4,",  heat. = ",e11.4)', &
   & Rate_energy_full_eVns1, &
   & Rate_energy_wall_eVns1(1) + Rate_energy_wall_eVns1(2), &
   & Rate_energy_emit_eVns1(1) + Rate_energy_emit_eVns1(2), &
   & Rate_energy_coll_eVns1(1) + Rate_energy_coll_eVns1(2), &
   & Rate_energy_heat_eVns1(1) + Rate_energy_heat_eVns1(2)

   PRINT  '("Sum of particular rates:     ",e11.4," eV/ns" )', &
   & Rate_energy_wall_eVns1(1) + Rate_energy_wall_eVns1(2) + &
   & Rate_energy_emit_eVns1(1) + Rate_energy_emit_eVns1(2) + &
   & Rate_energy_coll_eVns1(1) + Rate_energy_coll_eVns1(2) + &
   & Rate_energy_heat_eVns1(1) + Rate_energy_heat_eVns1(2)

   PRINT &
   & '(/,"Electrons, primary : number [1/ns] , <energy> [eV] : left = ",  e11.4," , ",  f11.4," ; right = ",  e11.4," , ",f11.4)',&
   & Rate_number_leftwall_ns1(1),       Avg_energy_leftwall_eV(1), &
   & Rate_number_rightwall_ns1(1),      Avg_energy_rightwall_eV(1)
   PRINT    '("                                  [macroparticles] : left = ",    i11,          14x," ; right = ", i11)', &
   & Rate_number_leftwall(1),           Rate_number_rightwall(1)
   PRINT &
   & '("Electrons, emitted : number [1/ns] , <energy> [eV] : left = ",  e11.4," , ",  f11.4," ; right = ",  e11.4," , ",f11.4)', &
   & Rate_number_leftemit_ns1(1),       Avg_energy_leftemit_eV(1), &
   & Rate_number_rightemit_ns1(1),      Avg_energy_rightemit_eV(1)
   PRINT &
   & '("                                  [macroparticles] : left = ",    i11,          14x," ; right = ", i11,13x," ::")', &
   & Rate_number_leftemit(1),           Rate_number_rightemit(1)

   see_right_elastic_energy   = GetAverageValue(see_right_elastic_energy,   see_right_elastic_count)   * Factor_energy_eV
   see_right_inelastic_energy = GetAverageValue(see_right_inelastic_energy, see_right_inelastic_count) * Factor_energy_eV
   see_right_true_energy      = GetAverageValue(see_right_true_energy,      see_right_true_count)      * Factor_energy_eV

   IF (PlasmaSourceFlag.EQ.0) THEN
      see_left_elastic_energy   = GetAverageValue(see_left_elastic_energy,   see_left_elastic_count)   * Factor_energy_eV
      see_left_inelastic_energy = GetAverageValue(see_left_inelastic_energy, see_left_inelastic_count) * Factor_energy_eV
      see_left_true_energy      = GetAverageValue(see_left_true_energy,      see_left_true_count)      * Factor_energy_eV

      !       SEE elastic : number [mcrp] , <energy> [eV] : left = 11111111111 , 11111111111 ; right = 11111111111 , 11111111111
      !     SEE inelastic : number [mcrp] , <energy> [eV] : left = 11111111111 , 11111111111 ; right = 11111111111 , 11111111111
      !          SEE true : number [mcrp] , <energy> [eV] : left = 11111111111 , 11111111111 ; right = 11111111111 , 11111111111
      !                                         SEE yield : left =     1111111               ; right =     1111111

      PRINT &
      & '("       SEE elastic : number [mcrp] , <energy> [eV] : left = ",    i11," , ",  f11.4," ; right = ",    i11," , ",f11.4)', &
      & see_left_elastic_count,        see_left_elastic_energy, &
      & see_right_elastic_count,       see_right_elastic_energy
      PRINT &
      & '("     SEE inelastic : number [mcrp] , <energy> [eV] : left = ",    i11," , ",  f11.4," ; right = ",    i11," , ",f11.4)', &
      & see_left_inelastic_count,      see_left_inelastic_energy, &
      & see_right_inelastic_count,     see_right_inelastic_energy
      PRINT &
      & '("          SEE true : number [mcrp] , <energy> [eV] : left = ",    i11," , ",  f11.4," ; right = ",    i11," , ",f11.4)', &
      & see_left_true_count,           see_left_true_energy, &
      & see_right_true_count,          see_right_true_energy
      PRINT '("                                         SEE yield : left =     ",f7.4,         14x," ; right =     ",f7.4)', &
      & Yield_left(1), Yield_right(1)
   ELSE
      electron_reflux_energy = GetAverageValue(electron_reflux_energy, electron_reflux_count) * Factor_energy_eV
      electron_reemit_energy = GetAverageValue(electron_reemit_energy, electron_reemit_count) * Factor_energy_eV
      !           refluxed, left : 11111111111 mcrp, 11111111111 eV ;   SEE elastic, right : 11111111111 mcrp, 11111111111 eV
      !         re-emitted, left : 11111111111 mcrp, 11111111111 eV ; SEE inelastic, right : 11111111111 mcrp, 11111111111 eV
      !                                                                    SEE true, right : 11111111111 mcrp, 11111111111 eV
      !    flux ratio (in/out), left = 1111111                      ;         SEE yield, right = 1111111

      PRINT &
      & '("           refluxed, left : ",    i11," mcrp, ",  f11.4," eV ;   SEE elastic, right : ",    i11," mcrp, ",  f11.4," eV")',&
      & electron_reflux_count,         electron_reflux_energy , &
      & see_right_elastic_count,       see_right_elastic_energy

      PRINT &
      & '("         re-emitted, left : ",    i11," mcrp, ",  f11.4," eV ; SEE inelastic, right : ",    i11," mcrp, ",  f11.4," eV")',&
      & electron_reemit_count,         electron_reemit_energy , &
      & see_right_inelastic_count,     see_right_inelastic_energy
      PRINT &
      & '("                                                                    SEE true, right : ",    i11," mcrp, ",  f11.4," eV")',&
      & see_right_true_count,          see_right_true_energy

      PRINT '("    flux ratio (in/out), left = ",f7.4,4x,"                  ;         SEE yield, right = ",f7.4)', &
      & Yield_left(1), Yield_right(1)
   END IF

   DO s = 2, N_spec
      !    Ions hit walls : number [1/ns] , <energy> [eV] : left = 11111111111 , 11111111111 ; right = 11111111111 , 11111111111
      !                                  [macroparticles] : left = 11111111111               ; right = 11111111111

      PRINT &
      & '(/,"    Ions hit walls : number [1/ns] , <energy> [eV] : left = ",  e11.4," , ",  f11.4," ; right = ",  e11.4," , ",f11.4)',&
      & Rate_number_leftwall_ns1(2),       Avg_energy_leftwall_eV(2), &
      & Rate_number_rightwall_ns1(2),      Avg_energy_rightwall_eV(2)
      PRINT    '("                                  [macroparticles] : left = ",    i11,          14x," ; right = ", i11)', &
      & Rate_number_leftwall(2),           Rate_number_rightwall(2)

      IF (PlasmaSourceFlag.NE.0) THEN
         ion_reflux_energy = GetAverageValue(ion_reflux_energy, ion_reflux_count) * MS(s) * Factor_energy_eV
         ion_reemit_energy = GetAverageValue(ion_reemit_energy, ion_reemit_count) * MS(s) * Factor_energy_eV

         !      Ions emitted : number [1/ns] , <energy> [eV] : left = 11111111111 , 11111111111 ; right = 11111111111 , 11111111111
         !                                  [macroparticles] : left = 11111111111               ; right = 11111111111               ::
         PRINT &
         & '("     Ions, emitted : number [1/ns] , <energy> [eV] : left = ",  e11.4," , ",  f11.4," ; right = ",  e11.4," , ",  f11.4)',&
         & Rate_number_leftemit_ns1(2),       Avg_energy_leftemit_eV(2), &
         & Rate_number_rightemit_ns1(2),      Avg_energy_rightemit_eV(2)
         PRINT &
         & '("                                  [macroparticles] : left = ",    i11,          14x," ; right = ", i11,13x," ::")', &
         & Rate_number_leftemit(2),           Rate_number_rightemit(2)
         !          refluxed : number [mcrp] , <energy> [eV] : left = 11111111111 , 11111111111
         !        re-emitted : number [mcrp] , <energy> [eV] : left = 11111111111 , 11111111111
         !                               flux ratio (in/out) : left =   111111111
         PRINT '("          refluxed : number [mcrp] , <energy> [eV] : left = ",    i11," , ",  f11.4)', &
         & ion_reflux_count, ion_reflux_energy
         PRINT '("        re-emitted : number [mcrp] , <energy> [eV] : left = ",    i11," , ",  f11.4)', &
         & ion_reemit_count, ion_reemit_energy
         PRINT '("                               flux ratio (in/out) : left =   ",f9.4)', Yield_left(2)

      ELSE IF (Ion_interac_model.EQ.1) THEN
         !      Ions emitted : number [1/ns] , <energy> [eV] : left = 11111111111 , 11111111111 ; right = 11111111111 , 11111111111
         !                                  [macroparticles] : left = 11111111111               ; right = 11111111111
         PRINT &
         & '("     Ions, emitted : number [1/ns] , <energy> [eV] : left = ",  e11.4," , ",  f11.4," ; right = ",  e11.4," , ",  f11.4)',&
         & Rate_number_leftemit_ns1(2),       Avg_energy_leftemit_eV(2), &
         & Rate_number_rightemit_ns1(2),      Avg_energy_rightemit_eV(2)
         PRINT '("                                  [macroparticles] : left = ",    i11,          14x," ; right = ", i11)', &
         & Rate_number_leftemit(2),           Rate_number_rightemit(2)

      END IF

   END DO

   PRINT '(/,"                  Average energy [eV]   :  el. = ",f11.4,",    ion = ",f11.4)', Avg_kin_energy_eV(1), &
   & Avg_kin_energy_eV(2)
   PRINT   '("                 Average density [m^-3] :  el. = ",e11.4,",    ion = ",e11.4)', N_avg_m3(1), N_avg_m3(2)
   PRINT   '("          Electron flow velocity [m/s]  : <VY> = ",f11.1,",   <VZ> = ",f11.1)', VY_e_avg_ms, VZ_e_avg_ms
   PRINT   '("Average electric current density [A/m^2]: <JY> = ",e11.4,",   <JZ> = ",e11.4)', JY_avg_Am2, JZ_avg_Am2

   IF (N_spec.EQ.1) THEN
      PRINT '(/"ELECTRONS, total macroparticles :",i9)', N_part(1)
   ELSE
      PRINT '(/"ELECTRONS, total macroparticles :",i9," ;  IONS, total macroparticles :",i9)', N_part(1), N_part(2)
   END IF

   PRINT &
   & '("e-n coll. : elastic-1 : ",i7,", <freq.> : ",e11.4," s^-1 ; i-n coll. : elastic-1: ",i7,", <freq.> : ",e11.4," s^-1")', &
   & e_n_1_count, e_n_1_count * f_factor(1), i_n_1_count, i_n_1_count * f_factor(2)
   PRINT &
   & '("         excitation-1 : ",i7,", <freq.> : ",e11.4," s^-1 ;     charge exchange-1: ",i7,", <freq.> : ",e11.4," s^-1")', &
   & e_n_2_count, e_n_2_count * f_factor(1), i_n_2_count, i_n_2_count * f_factor(2)
   PRINT &
   & '("         excitation-2 : ",i7,", <freq.> : ",e11.4," s^-1 ;")', e_n_5_count, e_n_5_count * f_factor(1)
   PRINT &
   & '("         ionization-1 : ",i7,", <freq.> : ",e11.4," s^-1 ;")', e_n_3_count, e_n_3_count * f_factor(1)
   PRINT &
   & '("         turbulence-1 : ",i7,", <freq.> : ",e11.4," s^-1 ;          turbulence-1: ",i7,", <freq.> : ",e11.4," s^-1")',&
   & e_t_4_count, e_t_4_count * f_factor(1), i_t_3_count, i_t_3_count * f_factor(2)

   balance_e = e_n_3_count + Rate_number_leftemit(1) + Rate_number_rightemit(1) - &
   & Rate_number_leftwall(1) - Rate_number_rightwall(1)
   balance_i = e_n_3_count + Rate_number_leftemit(2) + Rate_number_rightemit(2) - &
   & Rate_number_leftwall(2) - Rate_number_rightwall(2)

   PRINT '(/,"balance, electrons:",i7,", ions: ",i7,", total charge (i-e): ",i7)',  balance_e, balance_i, balance_i - balance_e

   PRINT '(//,"---------------------------------------------------------------------------------------------------------------")'

END SUBROUTINE DoTextOutput

!--------------------------------------
REAL(8) FUNCTION GetAverageValue(v, n)
!! Calculates average of double v and integer n as v/n, if n>0, and 0 otherwise.
!! Returns a double

   REAL(8) v
   INTEGER n

   if(n.EQ.0) then
      GetAverageValue = 0.0_8
   else
      GetAverageValue = v / n
   end if

END FUNCTION GetAverageValue

!--------------------------------------
REAL(8) FUNCTION GetAverageValueInt(m, n)
!! Calculates average for integers v and n as v/n, if n>0, and 0 otherwise.
!! Returns a double

   INTEGER m
   INTEGER n

   if(n.EQ.0) then
      GetAverageValueInt = 0.0_8
   else
      GetAverageValueInt = DBLE(m) / DBLE(n)
   end if

END FUNCTION GetAverageValueInt

!--------------------------------------------------
REAL(8) FUNCTION GetAverageFromBufer(length, bufer)

   IMPLICIT NONE

   INTEGER length
   REAL(8) bufer(1:length)
   INTEGER n, start, finish

   start = 1
   finish = length

   DO n = 1, length-2
! find first extremum
      IF ((bufer(n+1)-bufer(n))*(bufer(n+2)-bufer(n+1)).LE.0.0_8) THEN
         start = n+1
         EXIT
      END IF
   END DO

! find the next extremum
   DO n = start+1, length-2
      IF ((bufer(n+1)-bufer(n))*(bufer(n+2)-bufer(n+1)).LE.0.0_8) THEN
         finish = n+1
         EXIT
      END IF
   END DO

   GetAverageFromBufer = 0.0_8
   DO n = start, finish
      GetAverageFromBufer = GetAverageFromBufer + bufer(n)
   END DO
   GetAverageFromBufer = GetAverageFromBufer / MAX(1,finish-start+1)

END FUNCTION GetAverageFromBufer

!------------------------------
!
SUBROUTINE FINISH_DIAGNOSTICS

   USE ParallelOperationValues
   IMPLICIT NONE

   CALL DeallocateDiagnosticArrays

   IF (Rank_of_process.NE.0) RETURN

   PRINT *, ' Job finished! Terminate the program...'

END SUBROUTINE FINISH_DIAGNOSTICS

!-------------------------------------
!
SUBROUTINE AllocateDiagnosticArrays

   USE Diagnostics
   USE CurrentProblemValues, ONLY : N_cells, N_spec

   IMPLICIT NONE
   INTEGER ALLOC_ERR

! diagnostic arrays, cell-species
   IF (.NOT.ALLOCATED(VX2_cell)) THEN
      ALLOCATE(VX2_cell(1:N_cells,1:N_spec), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in ALLOCATE VX2_cell !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (.NOT.ALLOCATED(VY2_cell)) THEN
      ALLOCATE(VY2_cell(1:N_cells,1:N_spec), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in ALLOCATE VY2_cell !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (.NOT.ALLOCATED(VZ2_cell)) THEN
      ALLOCATE(VZ2_cell(1:N_cells,1:N_spec), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in ALLOCATE VZ2_cell !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (.NOT.ALLOCATED(VX_cell)) THEN
      ALLOCATE(VX_cell(1:N_cells,1:N_spec), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in ALLOCATE VX_cell !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (.NOT.ALLOCATED(VY_cell)) THEN
      ALLOCATE(VY_cell(1:N_cells,1:N_spec), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in ALLOCATE VY_cell !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (.NOT.ALLOCATED(VZ_cell)) THEN
      ALLOCATE(VZ_cell(1:N_cells,1:N_spec), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in ALLOCATE VZ_cell !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (.NOT.ALLOCATED(Npart_cell)) THEN
      ALLOCATE(Npart_cell(1:N_cells,1:N_spec), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in ALLOCATE Npart_cell !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

! diagnostic arrays, mesh
   IF (.NOT.ALLOCATED(QF_mesh)) THEN
      ALLOCATE(QF_mesh(0:N_cells), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in ALLOCATE QF_mesh !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (.NOT.ALLOCATED(QVX_mesh)) THEN
      ALLOCATE(QVX_mesh(0:N_cells), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in ALLOCATE QVX_mesh !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (.NOT.ALLOCATED(QVY_mesh)) THEN
      ALLOCATE(QVY_mesh(0:N_cells), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in ALLOCATE QVY_mesh !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (.NOT.ALLOCATED(QVZ_mesh)) THEN
      ALLOCATE(QVZ_mesh(0:N_cells), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in ALLOCATE QVZ_mesh !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (.NOT.ALLOCATED(NVX_mesh)) THEN
      ALLOCATE(NVX_mesh(0:N_cells, 1:N_spec), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in ALLOCATE NVX_mesh !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (.NOT. ALLOCATED(N_new_cell)) THEN
      ALLOCATE(N_new_cell(1:N_cells), STAT=ALLOC_ERR)
   END IF

   IF (.NOT. ALLOCATED(P_heat_cell)) THEN
      ALLOCATE(P_heat_cell(1:N_cells, 1:N_spec), STAT=ALLOC_ERR)
   END IF

END SUBROUTINE AllocateDiagnosticArrays

!---------------------------------------
!
SUBROUTINE DeallocateDiagnosticArrays

   USE Diagnostics
   IMPLICIT NONE
   INTEGER DEALLOC_ERR

! diagnostic arrays, cell-species
   IF (ALLOCATED(VX2_cell)) THEN
      DEALLOCATE(VX2_cell, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE VX2_cell !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(VY2_cell)) THEN
      DEALLOCATE(VY2_cell, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in ALLOCATE VY2_cell !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(VZ2_cell)) THEN
      DEALLOCATE(VZ2_cell, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE VZ2_cell !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(VX_cell)) THEN
      DEALLOCATE(VX_cell, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE VX_cell !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(VY_cell)) THEN
      DEALLOCATE(VY_cell, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE VY_cell !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(VZ_cell)) THEN
      DEALLOCATE(VZ_cell, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE VZ_cell !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(Npart_cell)) THEN
      DEALLOCATE(Npart_cell, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE Npart_cell !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

! diagnostic arrays, mesh
   IF (ALLOCATED(QF_mesh)) THEN
      DEALLOCATE(QF_mesh, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE QF_mesh !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(QVX_mesh)) THEN
      DEALLOCATE(QVX_mesh, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE QVX_mesh !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(QVY_mesh)) THEN
      DEALLOCATE(QVY_mesh, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE QVY_mesh !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(QVZ_mesh)) THEN
      DEALLOCATE(QVZ_mesh, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE QVZ_mesh !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(NVX_mesh)) THEN
      DEALLOCATE(NVX_mesh, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE NVX_mesh !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(N_new_cell)) THEN
      DEALLOCATE(N_new_cell, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE N_new_cell !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(P_heat_cell)) THEN
      DEALLOCATE(P_heat_cell, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE P_heat_cell !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

END SUBROUTINE DeallocateDiagnosticArrays
