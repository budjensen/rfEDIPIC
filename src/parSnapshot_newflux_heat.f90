
!------------------------------------------------------------------------------------------------------------
SUBROUTINE REQUEST_WALL_DF_FOR_SNAPSHOT

   USE Snapshots
   USE CurrentProblemValues, ONLY : T_cntr
   USE Diagnostics, ONLY : WriteAvg_step, WriteOut_step
   IMPLICIT NONE

!  INTEGER ALLOC_ERR

! quit if all snapshots were created or if due to any reasons counter of snapshots is larger than the declared number of snapshots (e.g. no snapshots)
   IF (current_snap.GT.N_of_all_snaps) RETURN

! if the current moment of time is the moment when it is necessary to start collecting velocities for velocity distribution function
! then set the corresponding flag value for .TRUE.

! Note, we turn on accumulation of particles in the timestep with number "Tcntr_snapshot(current_snap) - WriteAvg_step",
! while the diagnostic procedure DO_DIAGNOSTICS collects data starting at "Tcntr_snapshot(current_snap) - WriteAvg_step + 1" !
! We do so in order to account for the particles, which interact with the walls during the PREDICTED PUSH from
! "Tcntr_snapshot(current_snap) - WriteAvg_step" to "Tcntr_snapshot(current_snap) - WriteAvg_step + 1"

!  IF (T_cntr.EQ.(Tcntr_snapshot(current_snap) - WriteAvg_step)) THEN
   IF (T_cntr.GE.(Tcntr_snapshot(current_snap) - WriteOut_step)) THEN   ! we increased the interval of collecting
      Accumulate_wall_df = .TRUE.                        ! this flag is later turned off in subroutine CREATE_SNAPSHOT
      ! the corresponding distribution function arrays will be cleaned at that moment
      ! and therefore will be ready for a new accumulation cycle
   END IF

END SUBROUTINE REQUEST_WALL_DF_FOR_SNAPSHOT

!------------------------------------------------------------------------------------------------------------
!
SUBROUTINE CREATE_SNAPSHOT

   use mpi
   USE ParallelOperationValues
   USE Snapshots
   USE CurrentProblemValues
   USE Diagnostics
   USE LangevinCollisions
   IMPLICIT NONE

!  INCLUDE 'mpif.h'

   CHARACTER(19) evolF_filename
   CHARACTER(17) all_filename
   CHARACTER(16) Lg_filename

   INTEGER snap_blnks

   REAL(8) Ne_m3, Ni_m3             ! These are multiplied by 2 at the first and last node since N_scl_m3 = N_plasma_m3 / N_of_particles_cell.
                                    ! To tally up the total number of particles correctly we need to ensure each node point counts
                                    ! the particles in one cell, except for the end nodes which count the particles in half a cell
                                    ! (thus the number of particles is correctly tallied)
   REAL(8) Jx_Am2, Jy_Am2, Jz_Am2   ! These are multiplied by 2 at the first and last node since N_scl_m3 = N_plasma_m3 / N_of_particles_cell.
                                    ! To tally up the total number of particles correctly we need to ensure each node point counts
                                    ! the particles in one cell, except for the end nodes which count the particles in half a cell
                                    ! (thus the number of particles is correctly tallied)
   REAL(8) Vex_ms, Vey_ms, Vez_ms
   REAL(8) Wex_eV, Wey_eV, Wez_eV
   REAL(8) Vix_ms, Viy_ms, Viz_ms
   REAL(8) Wix_eV, Wiy_eV, Wiz_eV
   INTEGER n, i
   INTEGER s

   REAL(8) aa, bb, wxtemp, wytemp, wztemp !for calcualting "random" energy requested by Tim, 01/02/14
   REAL(8) density_e, ionization_sum, qe_m2s, qi_m2s
   REAL(8) factor_flux              !! factor_flux = coll_rate_factor / delta_t_s * N_in_macro
   REAL(8) factor_ionrate_int       !! factor_ionrate_int = factor_ionrate * delta_x_m
   REAL(8) factor_ionrate           !! factor_ionrate = (N_plasma_m3 / N_of_particles_cell) * coll_rate_factor / delta_t_s
   REAL(8) coll_rate_factor         !! Inverse of the collection time for flux counted by NVX_mesh
                                    !! coll_rate_factor = 1 / WriteOut_step
                                    !! OR               = Averaging_factor = 1 / WriteAvg_step for the first snapshot
   REAL(8) ionization_rate
   REAL(8) factor_Watt_cm3          !! Converts energy deposition due to collisions into physical units
                                    !! factor_Watt_cm3 = factor_ionrate * Factor_energy_eV * e_Cl * 1.e-6
   REAL(8) factor_eV                !! Unused conversion factor
                                    !! factor_eV = 0.5_8 * m_e_kg / e_Cl

   INTEGER ierr, ALLOC_ERR, DEALLOC_ERR
   INTEGER, ALLOCATABLE :: ibufer(:) !(1:N_cells)
   INTEGER, ALLOCATABLE :: ibufer2(:)
   REAL(8), ALLOCATABLE :: rbufer(:) !(1:N_nodes)
   REAL(8), ALLOCATABLE :: rbufer2(:)
   logical, parameter:: write_Lg = .false.
! functions
   REAL(8) Bx_gauss
   REAL(8) By_gauss

! check whether it is still necessary to save potential profiles
! note, counter_of_profiles_to_save can be non-zero only in the server process
   IF (counter_of_profiles_to_save.GT.0) THEN

! if we are here, then
! (i) this is the server node,
! (ii) the snapshot has been recently created and "current_snap" points to the next snapshot
!                                                (compared to the file number of device #90),
! but device #90 is still open for the server node

! so, we save the current time and the potential profile
      WRITE (90, '(1x,i8,1x,f11.5)') T_cntr, T_cntr * delta_t_s * 1.0d9
      WRITE (90, '(20(1x,f8.3))') F(0:N_cells) * U_scl_V

! update the profile counter
      counter_of_profiles_to_save = counter_of_profiles_to_save - 1

! close device #90 if all profiles were saved
      IF (counter_of_profiles_to_save.EQ.0) CLOSE (90, status = 'keep')

! and quit
      RETURN

! note, after counter_of_profiles_to_save becomes zero, we will not return here unless another snapshot will be created
   END IF

! quit if all snapshots were created or if due to any reasons counter of snapshots is larger than
! the declared number of snapshots (e.g. no snapshots)
   IF (current_snap.GT.N_of_all_snaps) RETURN

! quit if the current moment of time is not the moment when it is necessary to create the snapshot
   IF (T_cntr.NE.Tcntr_snapshot(current_snap)) RETURN

!###print '(" T_cntr = ",i8," process # ",i3," entered CREATE_SNAPSHOT and stayed within")', T_cntr, Rank_of_process

   IF (.NOT.ALLOCATED(ibufer)) THEN
      ALLOCATE(ibufer(1:N_cells), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : CREATE_SNAPSHOT : Error in ALLOCATE ibufer !!!")', Rank_of_process
         PRINT '(2x,"The program will be terminated now :(")'
         STOP
      END IF
   END IF

   IF (.NOT.ALLOCATED(ibufer2)) THEN
      ALLOCATE(ibufer2(1:N_cells), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : CREATE_SNAPSHOT : Error in ALLOCATE ibufer2 !!!")', Rank_of_process
         PRINT '(2x,"The program will be terminated now :(")'
         STOP
      END IF
   END IF

   IF (.NOT.ALLOCATED(rbufer)) THEN
      ALLOCATE(rbufer(1:N_nodes), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : CREATE_SNAPSHOT : Error in ALLOCATE rbufer !!!")', Rank_of_process
         PRINT '(2x,"The program will be terminated now :(")'
         STOP
      END IF
   END IF

   IF (.NOT.ALLOCATED(rbufer2)) THEN
      ALLOCATE(rbufer2(1:N_nodes), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : CREATE_SNAPSHOT : Error in ALLOCATE rbufer2 !!!")', Rank_of_process
         PRINT '(2x,"The program will be terminated now :(")'
         STOP
      END IF
   END IF

   IF (Rank_of_process.GT.0) THEN                                                                 ! client >>>
! send accumulated diagnostics to the server process
! N_new_cell
      DO s = 1, N_spec
         ibufer(1:N_cells)  = N_new_cell(1:N_cells)
         ibufer2(1:N_cells) = 0
         CALL MPI_REDUCE(ibufer, ibufer2, N_cells, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO
! Npart_cell
      DO s = 1, N_spec
         ibufer(1:N_cells)  = Npart_cell(1:N_cells, s)
         ibufer2(1:N_cells) = 0
         CALL MPI_REDUCE(ibufer, ibufer2, N_cells, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO
! P_heat_cell
      DO s = 1, N_spec
         rbufer(1:N_cells)  = P_heat_cell(1:N_cells, s)
         rbufer2(1:N_cells) = 0.0_8
         CALL MPI_REDUCE(rbufer, rbufer2, N_cells, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO
! VX2_cell
      DO s = 1, N_spec
         rbufer(1:N_cells)  = VX2_cell(1:N_cells, s)
         rbufer2(1:N_cells) = 0.0_8
         CALL MPI_REDUCE(rbufer, rbufer2, N_cells, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO
! VY2_cell
      DO s = 1, N_spec
         rbufer(1:N_cells)  = VY2_cell(1:N_cells, s)
         rbufer2(1:N_cells) = 0.0_8
         CALL MPI_REDUCE(rbufer, rbufer2, N_cells, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO
! VZ2_cell
      DO s = 1, N_spec
         rbufer(1:N_cells)  = VZ2_cell(1:N_cells, s)
         rbufer2(1:N_cells) = 0.0_8
         CALL MPI_REDUCE(rbufer, rbufer2, N_cells, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO
! VX_cell
      DO s = 1, N_spec
         rbufer(1:N_cells)  = VX_cell(1:N_cells, s)
         rbufer2(1:N_cells) = 0.0_8
         CALL MPI_REDUCE(rbufer, rbufer2, N_cells, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO
! VY_cell
      DO s = 1, N_spec
         rbufer(1:N_cells)  = VY_cell(1:N_cells, s)
         rbufer2(1:N_cells) = 0.0_8
         CALL MPI_REDUCE(rbufer, rbufer2, N_cells, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO
! VZ_cell
      DO s = 1, N_spec
         rbufer(1:N_cells)  = VZ_cell(1:N_cells, s)
         rbufer2(1:N_cells) = 0.0_8
         CALL MPI_REDUCE(rbufer, rbufer2, N_cells, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO
! QVX_mesh
      DO s = 1, N_spec
         rbufer(1:N_nodes)  = QVX_mesh(0:N_cells)
         rbufer2(1:N_nodes) = 0.0_8
         CALL MPI_REDUCE(rbufer, rbufer2, N_nodes, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO
! QVY_mesh
      DO s = 1, N_spec
         rbufer(1:N_nodes)  = QVY_mesh(0:N_cells)
         rbufer2(1:N_nodes) = 0.0_8
         CALL MPI_REDUCE(rbufer, rbufer2, N_nodes, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO
! QVZ_mesh
      DO s = 1, N_spec
         rbufer(1:N_nodes)  = QVZ_mesh(0:N_cells)
         rbufer2(1:N_nodes) = 0.0_8
         CALL MPI_REDUCE(rbufer, rbufer2, N_nodes, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO

! NVX_mesh
      DO s = 1, N_spec
         rbufer(1:N_nodes)  = NVX_mesh(0:N_cells,s)
         rbufer2(1:N_nodes) = 0.0_8
         CALL MPI_REDUCE(rbufer, rbufer2, N_nodes, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO

   ELSE
! receive data from clients

      PRINT '(/2x,"Process ",i3," : ^^^^^^^^^^^^^^^^^^^^ Snapshot ",i4," will be created now ... ^^^^^^^^^^^^^^^^^")', &
      & Rank_of_process, current_snap

! N_new_cell
      DO s = 1, N_spec
         ibufer  = 0
         ibufer2 = 0
         CALL MPI_REDUCE(ibufer2, ibufer, N_cells, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         N_new_cell(1:N_cells) = ibufer(1:N_cells)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO
! Npart_cell
      DO s = 1, N_spec
         ibufer  = 0
         ibufer2 = 0
         CALL MPI_REDUCE(ibufer2, ibufer, N_cells, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         Npart_cell(1:N_cells, s) = ibufer(1:N_cells)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO
! P_heat_cell
      DO s = 1, N_spec
         rbufer  = 0.0_8
         rbufer2 = 0.0_8
         CALL MPI_REDUCE(rbufer2, rbufer, N_cells, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         P_heat_cell(1:N_cells, s) = rbufer(1:N_cells)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO
! VX2_cell
      DO s = 1, N_spec
         rbufer  = 0.0_8
         rbufer2 = 0.0_8
         CALL MPI_REDUCE(rbufer2, rbufer, N_cells, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         VX2_cell(1:N_cells, s) = rbufer(1:N_cells)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO
! VY2_cell
      DO s = 1, N_spec
         rbufer  = 0.0_8
         rbufer2 = 0.0_8
         CALL MPI_REDUCE(rbufer2, rbufer, N_cells, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         VY2_cell(1:N_cells, s) = rbufer(1:N_cells)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO
! VZ2_cell
      DO s = 1, N_spec
         rbufer  = 0.0_8
         rbufer2 = 0.0_8
         CALL MPI_REDUCE(rbufer2, rbufer, N_cells, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         VZ2_cell(1:N_cells, s) = rbufer(1:N_cells)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO
! VX_cell
      DO s = 1, N_spec
         rbufer  = 0.0_8
         rbufer2 = 0.0_8
         CALL MPI_REDUCE(rbufer2, rbufer, N_cells, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         VX_cell(1:N_cells, s) = rbufer(1:N_cells)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO
! VY_cell
      DO s = 1, N_spec
         rbufer  = 0.0_8
         rbufer2 = 0.0_8
         CALL MPI_REDUCE(rbufer2, rbufer, N_cells, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         VY_cell(1:N_cells, s) = rbufer(1:N_cells)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO
! VZ_cell
      DO s = 1, N_spec
         rbufer  = 0.0_8
         rbufer2 = 0.0_8
         CALL MPI_REDUCE(rbufer2, rbufer, N_cells, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         VZ_cell(1:N_cells, s) = rbufer(1:N_cells)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO
! QVX_mesh
      DO s = 1, N_spec
         rbufer  = 0.0_8
         rbufer2 = 0.0_8
         CALL MPI_REDUCE(rbufer2, rbufer, N_nodes, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         QVX_mesh(0:N_cells) = rbufer(1:N_nodes)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO
! QVY_mesh
      DO s = 1, N_spec
         rbufer  = 0.0_8
         rbufer2 = 0.0_8
         CALL MPI_REDUCE(rbufer2, rbufer, N_nodes, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         QVY_mesh(0:N_cells) = rbufer(1:N_nodes)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO
! QVZ_mesh
      DO s = 1, N_spec
         rbufer  = 0.0_8
         rbufer2 = 0.0_8
         CALL MPI_REDUCE(rbufer2, rbufer, N_nodes, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         QVZ_mesh(0:N_cells) = rbufer(1:N_nodes)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO
! NVX_mesh
      DO s = 1, N_spec
         rbufer  = 0.0_8
         rbufer2 = 0.0_8
         CALL MPI_REDUCE(rbufer2, rbufer, N_nodes, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         NVX_mesh(0:N_cells,s) = rbufer(1:N_nodes)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO

! Produce new filenames
      WRITE (snapnumber_txt, '(i4)') current_snap       ! Conversation of integer i to character*4 number_txt
      snapnumber_txt = ADJUSTL(TRIM(snapnumber_txt))    ! Align to the left left in order to
      snap_blnks = 4 - LEN_TRIM(snapnumber_txt)         ! calculate the number of blanks in string number_txt;
      snapnumber_txt = ADJUSTR(snapnumber_txt)          ! then align to the right and fulfill
      snapnumber_txt(1:snap_blnks) = '0000'             ! substitution of leading blanks by '0'

! create the filenames

      evolF_filename =     '_TTTT_evolF_vsx.dat'
      evolF_filename(2:5) = snapnumber_txt

      all_filename =     '_TTTT_all_vsx.dat'
      all_filename(2:5) = snapnumber_txt

      Lg_filename =    '_TTTT_Lg_vsx.dat'
      Lg_filename(2:5) = snapnumber_txt

! write the data

!! evolution of the potential --------------------------------------------------------------------------------
!     OPEN (90, FILE = evolF_filename)
!! calculate the duration of the period of saving the potential profile's evolution
!     counter_of_profiles_to_save = 2.2_8 * &
!                                 & 2.0_8 * 3.1415926_8 * &
!                                 & SQRT( eps_0_Fm * m_e_kg * N_of_particles_cell / &
!                                 &       (Q_strm_spec(N_cells/2, 1) * N_plasma_m3 * e_Cl**2) ) &
!                                 & / delta_t_s
!
!!     print '(" #### uncorrected counter_of_profiles_to_save = ",i8)', counter_of_profiles_to_save
!
!     IF (current_snap.LT.N_of_all_snaps) THEN
!! cannot overlap with the next snapshot or be more than 1000
!        counter_of_profiles_to_save = MIN(counter_of_profiles_to_save, &
!                                         & 1000, &
!                                         & Tcntr_snapshot(current_snap+1)-Tcntr_snapshot(current_snap))
!     ELSE
!! cannot extend beyond the simulation time limit or be more than 1000
!        counter_of_profiles_to_save = MIN(counter_of_profiles_to_save, &
!                                         & 1000, &
!                                         & Max_T_cntr - T_cntr + 1)
!     END IF
!
!!     print '(" ####   corrected counter_of_profiles_to_save = ",i8)', counter_of_profiles_to_save
!
!! save the number of profiles (expected) and the number of points (assuming that point numbering starts from zero)
!     WRITE (90, '(1x,i4,1x,i4)') counter_of_profiles_to_save, N_cells
!! save the positions
!     WRITE (90, '(20(1x,f10.8))') (n*delta_x_m,n=0,N_cells,1)
!! save the current time
!     WRITE (90, '(1x,i8,1x,f11.5)') T_cntr, T_cntr * delta_t_s * 1.0d9
!! save the potential profile
!     WRITE (90, '(20(1x,f8.3))') F(0:N_cells) * U_scl_V
!! update the profile counter
!     counter_of_profiles_to_save = counter_of_profiles_to_save - 1
!! close device #90 if all profiles were saved
!     IF (counter_of_profiles_to_save.EQ.0) CLOSE (90, status = 'keep')
!! note, device #90 may remain open

!-------------------------------------------------------------------------------------------------------------
! note, in this file the values which previously were saved in half-integer nodes (middles of cells)
! are saved in the integer nodes (cells' boundaries)
      OPEN (99, FILE = all_filename)

! save column description
      WRITE (99, '("# col  1 is X-coordinate [m]")')

      WRITE (99, '("# col  2 is electrostatic potential [V]")')
      WRITE (99, '("# col  3 is X-electric field [V]")')

      WRITE (99, '("# col  4 is electron number density [m^-3]")')
      WRITE (99, '("# col  5 is ion number density [m^-3]")')

      WRITE (99, '("# col  6 is X-electric current density [A/m^2]")')
      WRITE (99, '("# col  7 is Y-electric current density [A/m^2]")')
      WRITE (99, '("# col  8 is Z-electric current density [A/m^2]")')

      WRITE (99, '("# col  9 is average X-velocity of electrons [m/s]")')
      WRITE (99, '("# col 10 is average Y-velocity of electrons [m/s]")')
      WRITE (99, '("# col 11 is average Z-velocity of electrons [m/s]")')

      WRITE (99, '("# col 12 is average X-energy of electrons [eV]")')
      WRITE (99, '("# col 13 is average Y-energy of electrons [eV]")')
      WRITE (99, '("# col 14 is average Z-energy of electrons [eV]")')

      WRITE (99, '("# col 15 is average X-velocity of ions [m/s]")')
      WRITE (99, '("# col 16 is average Y-velocity of ions [m/s]")')
      WRITE (99, '("# col 17 is average Z-velocity of ions [m/s]")')

      WRITE (99, '("# col 18 is average X-energy of ions [eV]")')
      WRITE (99, '("# col 19 is average Y-energy of ions [eV]")')
      WRITE (99, '("# col 20 is average Z-energy of ions [eV]")')

      WRITE (99, '("# col 21 is BX-magnetic field [Gauss]")')
      WRITE (99, '("# col 22 is BY-magnetic field [Gauss]")')

      WRITE (99, '("# col 23 is volume ionization rate Z(x)")')
      WRITE (99, '("# col 24 is integral of Z(x)dx)")')
      WRITE (99, '("# col 25 is ion flux")')
      WRITE (99, '("# col 26 is electron flux")')

      WRITE (99, '("# col 27 is gas heating rate in W/cm3 due to elastic collisions w. electrons")')
      WRITE (99, '("# col 28 is gas heating rate in W/cm3 due to collisions w. ions")')
      write (99, '("# col 29 is normalized charge density at x=0")')



! the values below are set zero for the case when there is only one species (N_spec=1)
      Ni_m3 = 0.0_8
      Vix_ms = 0.0_8
      Viy_ms = 0.0_8
      Viz_ms = 0.0_8
      Wix_eV = 0.0_8
      Wiy_eV = 0.0_8
      Wiz_eV = 0.0_8

      ionization_sum = 0.0_8
!!     factor_flux = Averaging_factor * (v_Te_ms / N_box_vel) * (N_plasma_m3 / N_of_particles_cell)
      factor_flux = Averaging_factor / delta_t_s * N_in_macro !flux calculated by counting particles crossing the node
      coll_rate_factor = 1. / dble(WriteOut_step)
      IF (current_snap .EQ. 1) coll_rate_factor = Averaging_factor
      factor_flux = coll_rate_factor / delta_t_s * N_in_macro
      factor_ionrate = (N_plasma_m3 / N_of_particles_cell) * coll_rate_factor / delta_t_s
      factor_ionrate_int = factor_ionrate * delta_x_m
      factor_Watt_cm3 = factor_ionrate * Factor_energy_eV * e_Cl * 1.e-6 !convert energy deposition due to collisions into phys. units

      DO n = 0, N_cells

         IF (n.EQ.0) THEN
            aa = MAX(DBLE(Npart_cell(1,1)),1.0_8)
            Ne_m3 = 2.0_8 * Q_strm_spec(n, 1) * N_scl_m3
            Jx_Am2 = 2.0_8 * QVX_mesh(n) * J_scl_Am2
            Jy_Am2 = 2.0_8 * QVY_mesh(n) * J_scl_Am2
            Jz_Am2 = 2.0_8 * QVZ_mesh(n) * J_scl_Am2
            Vex_ms = (VX_cell(1,1) / MAX(DBLE(Npart_cell(1,1)),1.0_8)) * V_scl_ms
            Vey_ms = (VY_cell(1,1) / MAX(DBLE(Npart_cell(1,1)),1.0_8)) * V_scl_ms
            Vez_ms = (VZ_cell(1,1) / MAX(DBLE(Npart_cell(1,1)),1.0_8)) * V_scl_ms
            Wex_eV = (VX2_cell(1,1) / MAX(DBLE(Npart_cell(1,1)),1.0_8)) * Factor_energy_eV * Ms(1)
            Wey_eV = (VY2_cell(1,1) / MAX(DBLE(Npart_cell(1,1)),1.0_8)) * Factor_energy_eV * Ms(1)
            Wez_eV = (VZ2_cell(1,1) / MAX(DBLE(Npart_cell(1,1)),1.0_8)) * Factor_energy_eV * Ms(1)
!           Wex_eV = (VX2_cell(1,1)/aa - (VX_cell(1,1) / aa)**2) * Factor_energy_eV * Ms(1)
!           Wey_eV = (VY2_cell(1,1)/aa - (VY_cell(1,1) / aa)**2) * Factor_energy_eV * Ms(1)
!           Wez_eV = (VZ2_cell(1,1)/aa - (VZ_cell(1,1) / aa)**2) * Factor_energy_eV * Ms(1)
            qe_m2s = NVX_mesh(n, 1) * factor_flux
            ionization_rate = 0.5 * N_new_cell(1) * factor_ionrate
            ionization_sum = 0.5 * N_new_cell(1) * factor_ionrate_int
            IF (N_spec.EQ.2) THEN
               aa = MAX(DBLE(Npart_cell(1,2)),1.0_8)
               Ni_m3 = 2.0_8 * Q_strm_spec(n, 2) * N_scl_m3
               Vix_ms = (VX_cell(1,2) / MAX(DBLE(Npart_cell(1,2)),1.0_8)) * V_scl_ms
               Viy_ms = (VY_cell(1,2) / MAX(DBLE(Npart_cell(1,2)),1.0_8)) * V_scl_ms
               Viz_ms = (VZ_cell(1,2) / MAX(DBLE(Npart_cell(1,2)),1.0_8)) * V_scl_ms
               Wix_eV = (VX2_cell(1,2) / MAX(DBLE(Npart_cell(1,2)),1.0_8)) * Factor_energy_eV * Ms(2)
               Wiy_eV = (VY2_cell(1,2) / MAX(DBLE(Npart_cell(1,2)),1.0_8)) * Factor_energy_eV * Ms(2)
               Wiz_eV = (VZ2_cell(1,2) / MAX(DBLE(Npart_cell(1,2)),1.0_8)) * Factor_energy_eV * Ms(2)
!              Wix_eV = (VX2_cell(1,2)/aa - (VX_cell(1,2) / aa)**2) * Factor_energy_eV * Ms(2)
!              Wiy_eV = (VY2_cell(1,2)/aa - (VY_cell(1,2) / aa)**2) * Factor_energy_eV * Ms(2)
!              Wiz_eV = (VZ2_cell(1,2)/aa - (VZ_cell(1,2) / aa)**2) * Factor_energy_eV * Ms(2)
               qi_m2s = NVX_mesh(n, 2) * factor_flux
            END IF
         ELSE IF (n.EQ.N_cells) THEN
            aa = MAX(DBLE(Npart_cell(n,1)),1.0_8)
            Ne_m3 = 2.0_8 * Q_strm_spec(n, 1) * N_scl_m3
            Jx_Am2 = 2.0_8 * QVX_mesh(n) * J_scl_Am2
            Jy_Am2 = 2.0_8 * QVY_mesh(n) * J_scl_Am2
            Jz_Am2 = 2.0_8 * QVZ_mesh(n) * J_scl_Am2
            Vex_ms = (VX_cell(n,1) / MAX(DBLE(Npart_cell(n,1)),1.0_8)) * V_scl_ms
            Vey_ms = (VY_cell(n,1) / MAX(DBLE(Npart_cell(n,1)),1.0_8)) * V_scl_ms
            Vez_ms = (VZ_cell(n,1) / MAX(DBLE(Npart_cell(n,1)),1.0_8)) * V_scl_ms
            Wex_eV = (VX2_cell(n,1) / MAX(DBLE(Npart_cell(n,1)),1.0_8)) * Factor_energy_eV * Ms(1)
            Wey_eV = (VY2_cell(n,1) / MAX(DBLE(Npart_cell(n,1)),1.0_8)) * Factor_energy_eV * Ms(1)
            Wez_eV = (VZ2_cell(n,1) / MAX(DBLE(Npart_cell(n,1)),1.0_8)) * Factor_energy_eV * Ms(1)
!           Wex_eV = (VX2_cell(n,1)/aa - (VX_cell(n,1) / aa)**2) * Factor_energy_eV * Ms(1)
!           Wey_eV = (VY2_cell(n,1)/aa - (VY_cell(n,1) / aa)**2) * Factor_energy_eV * Ms(1)
!           Wez_eV = (VZ2_cell(n,1)/aa - (VZ_cell(n,1) / aa)**2) * Factor_energy_eV * Ms(1)
            qe_m2s = NVX_mesh(n, 1) * factor_flux
            ionization_rate = 0.5 * N_new_cell(n) * factor_ionrate
            ionization_sum = ionization_sum + 0.5 * N_new_cell(n) * factor_ionrate_int
            IF (N_spec.EQ.2) THEN
               aa = MAX(DBLE(Npart_cell(n,2)),1.0_8)
               Ni_m3 = 2.0_8 * Q_strm_spec(n, 2) * N_scl_m3
               Vix_ms = (VX_cell(n,2) / MAX(DBLE(Npart_cell(n,2)),1.0_8)) * V_scl_ms
               Viy_ms = (VY_cell(n,2) / MAX(DBLE(Npart_cell(n,2)),1.0_8)) * V_scl_ms
               Viz_ms = (VZ_cell(n,2) / MAX(DBLE(Npart_cell(n,2)),1.0_8)) * V_scl_ms
               Wix_eV = (VX2_cell(n,2) / MAX(DBLE(Npart_cell(n,2)),1.0_8)) * Factor_energy_eV * Ms(2)
               Wiy_eV = (VY2_cell(n,2) / MAX(DBLE(Npart_cell(n,2)),1.0_8)) * Factor_energy_eV * Ms(2)
               Wiz_eV = (VZ2_cell(n,2) / MAX(DBLE(Npart_cell(n,2)),1.0_8)) * Factor_energy_eV * Ms(2)
!              Wix_eV = (VX2_cell(n,2)/aa - (VX_cell(n,2) / aa)**2) * Factor_energy_eV * Ms(2)
!              Wiy_eV = (VY2_cell(n,2)/aa - (VY_cell(n,2) / aa)**2) * Factor_energy_eV * Ms(2)
!              Wiz_eV = (VZ2_cell(n,2)/aa - (VZ_cell(n,2) / aa)**2) * Factor_energy_eV * Ms(2)
               qi_m2s = NVX_mesh(n, 2) * factor_flux
            END IF
         ELSE
            aa = MAX(DBLE(Npart_cell(n,1)),1.0_8)
            bb = MAX(DBLE(Npart_cell(n+1,1)),1.0_8)
            wxtemp = 0.5_8 * ( VX2_cell(n,1)/aa + VX2_cell(n+1,1)/bb )
            wytemp = 0.5_8 * ( VY2_cell(n,1)/aa + VY2_cell(n+1,1)/bb )
            wztemp = 0.5_8 * ( VZ2_cell(n,1)/aa + VZ2_cell(n+1,1)/bb )
            Ne_m3 = Q_strm_spec(n, 1) * N_scl_m3
            Jx_Am2 = QVX_mesh(n) * J_scl_Am2
            Jy_Am2 = QVY_mesh(n) * J_scl_Am2
            Jz_Am2 = QVZ_mesh(n) * J_scl_Am2
            Vex_ms = 0.5_8 * ( VX_cell(n,1)/MAX(DBLE(Npart_cell(n,1)),1.0_8) + VX_cell(n+1,1)/MAX(DBLE(Npart_cell(n+1,1)),1.0_8) ) * V_scl_ms
            Vey_ms = 0.5_8 * ( VY_cell(n,1)/MAX(DBLE(Npart_cell(n,1)),1.0_8) + VY_cell(n+1,1)/MAX(DBLE(Npart_cell(n+1,1)),1.0_8) ) * V_scl_ms
            Vez_ms = 0.5_8 * ( VZ_cell(n,1)/MAX(DBLE(Npart_cell(n,1)),1.0_8) + VZ_cell(n+1,1)/MAX(DBLE(Npart_cell(n+1,1)),1.0_8) ) * V_scl_ms
            Wex_eV = 0.5_8 * ( VX2_cell(n,1)/MAX(DBLE(Npart_cell(n,1)),1.0_8) + VX2_cell(n+1,1)/MAX(DBLE(Npart_cell(n+1,1)),1.0_8) ) * Factor_energy_eV * Ms(1)
            Wey_eV = 0.5_8 * ( VY2_cell(n,1)/MAX(DBLE(Npart_cell(n,1)),1.0_8) + VY2_cell(n+1,1)/MAX(DBLE(Npart_cell(n+1,1)),1.0_8) ) * Factor_energy_eV * Ms(1)
            Wez_eV = 0.5_8 * ( VZ2_cell(n,1)/MAX(DBLE(Npart_cell(n,1)),1.0_8) + VZ2_cell(n+1,1)/MAX(DBLE(Npart_cell(n+1,1)),1.0_8) ) * Factor_energy_eV * Ms(1)
!           Wex_eV = 0.5_8 * ( 2.0_8 * wxtemp - (VX_cell(n,1)/aa)**2 - (VX_cell(n+1,1)/bb)**2 ) * Factor_energy_eV * Ms(1)
!           Wey_eV = 0.5_8 * ( 2.0_8 * wytemp - (VY_cell(n,1)/aa)**2 - (VY_cell(n+1,1)/bb)**2 ) * Factor_energy_eV * Ms(1)
!           Wez_eV = 0.5_8 * ( 2.0_8 * wztemp - (VZ_cell(n,1)/aa)**2 - (VZ_cell(n+1,1)/bb)**2 ) * Factor_energy_eV * Ms(1)
            qe_m2s = NVX_mesh(n, 1) * factor_flux
            ionization_rate = 0.5 * (N_new_cell(n)+N_new_cell(n+1)) * factor_ionrate
            ionization_sum = ionization_sum + 0.5 * (N_new_cell(n) + N_new_cell(n+1)) * factor_ionrate_int
            IF (N_spec.EQ.2) THEN
               aa = MAX(DBLE(Npart_cell(n,2)),1.0_8)
               bb = MAX(DBLE(Npart_cell(n+1,2)),1.0_8)
               wxtemp = 0.5_8 * ( VX2_cell(n,2)/aa + VX2_cell(n+1,2)/bb )
               wytemp = 0.5_8 * ( VY2_cell(n,2)/aa + VY2_cell(n+1,2)/bb )
               wztemp = 0.5_8 * ( VZ2_cell(n,2)/aa + VZ2_cell(n+1,2)/bb )
               Ni_m3 = Q_strm_spec(n, 2) * N_scl_m3
               Vix_ms = 0.5_8 * ( VX_cell(n,2)/MAX(DBLE(Npart_cell(n,2)),1.0_8) + VX_cell(n+1,2)/MAX(DBLE(Npart_cell(n+1,2)),1.0_8) ) * V_scl_ms
               Viy_ms = 0.5_8 * ( VY_cell(n,2)/MAX(DBLE(Npart_cell(n,2)),1.0_8) + VY_cell(n+1,2)/MAX(DBLE(Npart_cell(n+1,2)),1.0_8) ) * V_scl_ms
               Viz_ms = 0.5_8 * ( VZ_cell(n,2)/MAX(DBLE(Npart_cell(n,2)),1.0_8) + VZ_cell(n+1,2)/MAX(DBLE(Npart_cell(n+1,2)),1.0_8) ) * V_scl_ms
               Wix_eV = 0.5_8 * ( VX2_cell(n,2)/MAX(DBLE(Npart_cell(n,2)),1.0_8) + VX2_cell(n+1,2)/MAX(DBLE(Npart_cell(n+1,2)),1.0_8) ) * Factor_energy_eV * Ms(2)
               Wiy_eV = 0.5_8 * ( VY2_cell(n,2)/MAX(DBLE(Npart_cell(n,2)),1.0_8) + VY2_cell(n+1,2)/MAX(DBLE(Npart_cell(n+1,2)),1.0_8) ) * Factor_energy_eV * Ms(2)
               Wiz_eV = 0.5_8 * ( VZ2_cell(n,2)/MAX(DBLE(Npart_cell(n,2)),1.0_8) + VZ2_cell(n+1,2)/MAX(DBLE(Npart_cell(n+1,2)),1.0_8) ) * Factor_energy_eV * Ms(2)
!              Wix_eV = 0.5_8 * ( 2.0_8 * wxtemp - (VX_cell(n,2)/aa)**2 - (VX_cell(n+1,2)/bb)**2 ) * Factor_energy_eV * Ms(2)
!              Wiy_eV = 0.5_8 * ( 2.0_8 * wytemp - (VY_cell(n,2)/aa)**2 - (VY_cell(n+1,2)/bb)**2 ) * Factor_energy_eV * Ms(2)
!              Wiz_eV = 0.5_8 * ( 2.0_8 * wztemp - (VZ_cell(n,2)/aa)**2 - (VZ_cell(n+1,2)/bb)**2 ) * Factor_energy_eV * Ms(2)
               qi_m2s = NVX_mesh(n, 2) * factor_flux
            END IF
         END IF

         factor_eV = 0.5_8 * m_e_kg / e_Cl

         IF (n.EQ.0) THEN
            WRITE (99, '(29(1x,e14.7))') &
            & n * delta_x_m, &   ! 1
            & F(n) * U_scl_V, &       ! 2
            & EX(n) * E_scl_Vm, &     ! 3
            & Ne_m3, &           ! 4
            & Ni_m3, &           ! 5
!              & Jx_Am2, &               ! 6
            & (qi_m2s - qe_m2s) * e_Cl, & ! 6
            & Jy_Am2, &               ! 7
            & Jz_Am2, &               ! 8
            & Vex_ms, &          !  9
            & Vey_ms, &          ! 10
            & Vez_ms, &          ! 11
!              & Wex_eV - Vex_ms**2 * Ms(1) * factor_eV, &     ! 12
!              & Wey_eV - Vey_ms**2 * Ms(1) * factor_eV, &     ! 13
!              & Wez_eV - Vez_ms**2 * Ms(1) * factor_eV, &     ! 14
            & Wex_eV, &      ! 12
            & Wey_eV, &      ! 13
            & Wez_eV, &      ! 14
            & Vix_ms, &          ! 15
            & Viy_ms, &          ! 16
            & Viz_ms, &          ! 17
!              & Wix_eV - Vix_ms**2 * Ms(2) * factor_eV, &     ! 18
!              & Wiy_eV - Viy_ms**2 * Ms(2) * factor_eV, &     ! 19
!              & Wiz_eV - Viz_ms**2 * Ms(2) * factor_eV, &     ! 20
            & Wix_eV, &     ! 18
            & Wiy_eV, &     ! 19
            & Wiz_eV, &     ! 20
            & Bx_gauss(DBLE(n)), &   ! 21
            & By_gauss(DBLE(n)), &   ! 22
            & ionization_rate, &     ! 23
            & ionization_sum, &      ! 24
            & qi_m2s, &              ! 25
            & qe_m2s, &              ! 26
            & 0.0_8, & !27 !!!! 9-22-23 - This use to print n=0, which was not allocated
            & 0.0_8, & !28 !!!! Now we just print P_heat_cell(n=0)=0, which is not a big deal (I hope)
            & full_Q_left                                    !29
         ELSE

            WRITE (99, '(29(1x,e14.7))') &
            & n * delta_x_m, &   ! 1
            & F(n) * U_scl_V, &       ! 2
            & EX(n) * E_scl_Vm, &     ! 3
            & Ne_m3, &           ! 4
            & Ni_m3, &           ! 5
            !             & Jx_Am2, &               ! 6
            & (qi_m2s - qe_m2s) * e_Cl, & ! 6
            & Jy_Am2, &               ! 7
            & Jz_Am2, &               ! 8
            & Vex_ms, &          !  9
            & Vey_ms, &          ! 10
            & Vez_ms, &          ! 11
            !             & Wex_eV - Vex_ms**2 * Ms(1) * factor_eV, &     ! 12
            !             & Wey_eV - Vey_ms**2 * Ms(1) * factor_eV, &     ! 13
            !             & Wez_eV - Vez_ms**2 * Ms(1) * factor_eV, &     ! 14
            & Wex_eV, &      ! 12
            & Wey_eV, &      ! 13
            & Wez_eV, &      ! 14
            & Vix_ms, &          ! 15
            & Viy_ms, &          ! 16
            & Viz_ms, &          ! 17
            !             & Wix_eV - Vix_ms**2 * Ms(2) * factor_eV, &     ! 18
            !             & Wiy_eV - Viy_ms**2 * Ms(2) * factor_eV, &     ! 19
            !             & Wiz_eV - Viz_ms**2 * Ms(2) * factor_eV, &     ! 20
            & Wix_eV, &     ! 18
            & Wiy_eV, &     ! 19
            & Wiz_eV, &     ! 20
            & Bx_gauss(DBLE(n)), &   ! 21
            & By_gauss(DBLE(n)), &   ! 22
            & ionization_rate, &     ! 23
            & ionization_sum, &      ! 24
            & qi_m2s, &              ! 25
            & qe_m2s, &              ! 26
            & P_heat_cell(n, 1) * Ms(1) * factor_Watt_cm3, & !27 !!!! 9-22-23 - This use to print n=0, which was not allocated
            & P_heat_cell(n, 2) * Ms(2) * factor_Watt_cm3, & !28 !!!! Now we just print P_heat_cell(n=0)=0, which is not a big deal (I hope)
            & full_Q_left                                    !29
         END IF

      END DO
      CLOSE (99, STATUS = 'KEEP')

! Langevin coefficients ------------------------------------------------------------------------------------------
! for e-e collisions only
      IF (Accounted_Langevin.EQ.1 .and. write_Lg) THEN
         OPEN (99, FILE = Lg_filename)
! save column description
         WRITE (99, '("# col 1 is X-coordinate [m]")')
         WRITE (99, '("# col 2 is absolute velocity [in units of scale electron thermal speed V_th_e = ",e12.5," m/s]")') V_Te_ms
         WRITE (99, '("#-----------------------------")')
         WRITE (99, '("# col 3 is the electron distribution function over the absolute velocity [arb.units]")')
         WRITE (99, '("# col 4 is the typical electron velocity change due to e-e DRAG FORCE in the direction")')
         WRITE (99, '("#          along the electron velocity in the electron flow frame [same units as col 2]")')
         WRITE (99, '("# col 5 is the typical electron velocity change ACROSS the electron velocity in the electron flow frame")')
         WRITE (99, '("#          due to the velocity DIFFUSION caused by e-e collisions [same units as col 2]")')
         WRITE (99, '("# col 6 is the typical electron velocity change ALONG the electron velocity in the electron flow frame")')
         WRITE (99, '("#          due to the velocity DIFFUSION caused by e-e collisions [same units as col 2]")')
         DO n = 0, N_cells-1

            IF (n.EQ.0) THEN
               density_e = 0.5_8 * (2.0_8 * Q_strm_spec(n,1) + Q_strm_spec(n+1,1))
            ELSE IF (n.EQ.N_cells-1) THEN
               density_e = 0.5_8 * (Q_strm_spec(n,1) + 2.0_8 * Q_strm_spec(n+1,1))
            ELSE
               density_e = 0.5_8 * (Q_strm_spec(n,1) + Q_strm_spec(n+1,1))
            END IF

            DO i = 0, N_box_w_Lang - 1
               WRITE (99, '(6(2x,e14.7))') &
               & n * delta_x_m, &                  ! 1
               & DBLE(i) / DBLE(N_box_vel), &            ! 2
               &     Fd_w(i+1, n), &                          ! 3
               &   F_drag(i,   n) * density_e / DBLE(N_box_vel), &               ! 4
               & D_1_sqrt(i,   n) * SQRT(density_e), &               ! 5 here we omitted N_box_vel/N_box_vel
               & D_3_sqrt(i,   n) * SQRT(density_e)                  ! 6 in numeriator, N_box_vel is from GetMaxwellVelocity
            END DO                                                        !   in denominator, N_box_vel is  from the velocity scaling

            WRITE (99, '(7(2x,e14.7))') &
            & n * delta_x_m, &                  ! 1
            & DBLE(N_box_w_Lang) / DBLE(N_box_vel), & ! 2
            & 0.0_8, &                                     ! 3
            &   F_drag(N_box_w_Lang, n) * density_e / DBLE(N_box_vel), &                 ! 4
            & D_1_sqrt(N_box_w_Lang, n) * SQRT(density_e), &                 ! 5
            & D_3_sqrt(N_box_w_Lang, n) * SQRT(density_e), &                    ! 6
            & Log_coul(n)                                                       ! 7

            WRITE (99, '(" ")')

         END DO
         CLOSE (99, STATUS = 'KEEP')
      END IF

! for e-e and e-i collisions
      IF (Accounted_Langevin.EQ.2 .and. write_Lg) THEN

         OPEN (99, FILE = Lg_filename)

! save column description
         WRITE (99, '("# col 1 is X-coordinate [m]")')
         WRITE (99, '("# col 2 is absolute velocity [in units of scale electron thermal speed V_th_e = ",e12.5," m/s]")') V_Te_ms
         WRITE (99, '("#-----------------------------")')
         WRITE (99, '("# col 3 is the electron distribution function over the absolute velocity [arb.units]")')
         WRITE (99, '("# col 4 is the typical electron velocity change due to e-e DRAG FORCE in the direction")')
         WRITE (99, '("#          along the electron velocity in the electron flow frame [same units as col 2]")')
         WRITE (99, '("# col 5 is the typical electron velocity change ACROSS the electron velocity in the electron flow frame")')
         WRITE (99, '("#          due to the velocity DIFFUSION caused by e-e collisions [same units as col 2]")')
         WRITE (99, '("# col 6 is the typical electron velocity change ALONG the electron velocity in the electron flow frame")')
         WRITE (99, '("#          due to the velocity DIFFUSION caused by e-e collisions [same units as col 2]")')
         WRITE (99, '("#-----------------------------")')
         WRITE (99, '("# col 7 is the typical electron velocity change due to e-i DRAG FORCE in the direction")')
         WRITE (99, '("#          along the electron velocity (the ion flow speed is neglected) [same units as col 2]")')
         WRITE (99, '("# col 8 is the typical electron velocity change ACROSS the electron velocity in the lab frame")')
         WRITE (99, '("#          due to the velocity DIFFUSION caused by e-i collisions [same units as col 2]")')
         DO n = 0, N_cells-1

            IF (n.EQ.0) THEN
               density_e = 0.5_8 * (2.0_8 * Q_strm_spec(n,1) + Q_strm_spec(n+1,1))
            ELSE IF (n.EQ.N_cells-1) THEN
               density_e = 0.5_8 * (Q_strm_spec(n,1) + 2.0_8 * Q_strm_spec(n+1,1))
            ELSE
               density_e = 0.5_8 * (Q_strm_spec(n,1) + Q_strm_spec(n+1,1))
            END IF

            DO i = 0, N_box_w_Lang - 1
               WRITE (99, '(8(2x,e14.7))') &
               & n * delta_x_m, &                  ! 1
               & DBLE(i) / DBLE(N_box_vel), &            ! 2
               &     Fd_w(i+1, n), &                          ! 3
               &   F_drag(i,   n) * density_e / DBLE(N_box_vel), &               ! 4
               & D_1_sqrt(i,   n) * SQRT(density_e), &               ! 5
               & D_3_sqrt(i,   n) * SQRT(density_e), &               ! 6
               & ( -0.039788736_8 * Factor_F_drag   *      Log_coul(n)  *      density_e  /     (MAX(DBLE(i), V_ion_threshold(n)))**2 ) / DBLE(N_box_vel), &  ! 7
               &    0.488602512_8 * Factor_D_1_sqrt * SQRT(Log_coul(n)) * SQRT(density_e) / SQRT(MAX(DBLE(i), V_ion_threshold(n)))        ! 8  similar to columns 5 and 6, we
            END DO                                                                                                                             !    omitted N_box_vel / N_box_vel

            WRITE (99, '(8(2x,e14.7))') &
            & n * delta_x_m, &                  ! 1
            & DBLE(N_box_w_Lang) / DBLE(N_box_vel), & ! 2
            & 0.0_8, &                                     ! 3
            &   F_drag(N_box_w_Lang, n) * density_e / DBLE(N_box_vel), &                 ! 4
            & D_1_sqrt(N_box_w_Lang, n) * SQRT(density_e), &                 ! 5
            & D_3_sqrt(N_box_w_Lang, n) * SQRT(density_e), &                 ! 6
            & ( -0.039788736_8 * Factor_F_drag   *      Log_coul(n)  *      density_e  /     (MAX(DBLE(N_box_w_Lang), V_ion_threshold(n)))**2 ) / DBLE(N_box_vel), &  ! 7
            &    0.488602512_8 * Factor_D_1_sqrt * SQRT(Log_coul(n)) * SQRT(density_e) / SQRT(MAX(DBLE(N_box_w_Lang), V_ion_threshold(n)))        ! 8

            WRITE (99, '(" ")')

         END DO

         CLOSE (99, STATUS = 'KEEP')

      END IF

   END IF

!!** disabled  CALL SNAP_ELECTRON_2D_VDF_WALLS
   if (flag_evxdf.or.flag_evydf.or.flag_evzdf) CALL SNAP_LOCAL_ELECTRON_VDFS
   if(flag_eedf) CALL SNAP_LOCAL_ELECTRON_EDFS
   CALL SNAP_ELECTRON_PHASE_PLANES

   IF (N_spec.EQ.2) THEN

      if (flag_ivxdf) CALL SNAP_LOCAL_ION_VDFS
      if (flag_iedf) CALL SNAP_LOCAL_ION_EDFS
      CALL SNAP_ION_PHASE_PLANE
      if (flag_ilwedf.or.flag_irwedf) CALL SNAP_ION_EDF_RW_LW

   END IF

   IF (Rank_of_process.EQ.0) PRINT &
   & '(/2x,"Process ",i3," : ^^^^^^^^^^^^^^^^^^^^ Snapshot ",i4," completed :) ^^^^^^^^^^^^^^^^^^^")', &
   & Rank_of_process, current_snap

   Accumulate_wall_df = .FALSE.              ! disable collecting of velocities for creation of wall distributions by SNAPSHOT procedures
   current_snap = current_snap + 1           ! increase the snapshots counter

   IF (ALLOCATED(ibufer)) THEN
      DEALLOCATE(ibufer, STAT = DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : CREATE_SNAPSHOT : Error in DEALLOCATE ibufer !!!")', Rank_of_process
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(ibufer2)) THEN
      DEALLOCATE(ibufer2, STAT = DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : CREATE_SNAPSHOT : Error in DEALLOCATE ibufer2 !!!")', Rank_of_process
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(rbufer)) THEN
      DEALLOCATE(rbufer, STAT = DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : CREATE_SNAPSHOT : Error in DEALLOCATE rbufer !!!")', Rank_of_process
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(rbufer2)) THEN
      DEALLOCATE(rbufer2, STAT = DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : CREATE_SNAPSHOT : Error in DEALLOCATE rbufer2 !!!")', Rank_of_process
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

END SUBROUTINE CREATE_SNAPSHOT

!-----------------------------------------------------------------------------------------------
! produces the 2-d velocity distribution functions for electrons hitting the right wall
!
SUBROUTINE SNAP_ELECTRON_2D_VDF_WALLS

   use mpi
   USE ParallelOperationValues
   USE Snapshots
   USE CurrentProblemValues !, ONLY : V_Te_ms

   IMPLICIT NONE

!  INCLUDE 'mpif.h'

   INTEGER j               ! index of a box of velocity parallel to the walls
   INTEGER i               ! index of a box of velocity normal to the walls

   CHARACTER(16) e2vdfw_filename

   INTEGER, ALLOCATABLE :: ibufer(:) !(0:N_box_Vx_e)
   INTEGER, ALLOCATABLE :: ibufer2(:)
   INTEGER ALLOC_ERR, DEALLOC_ERR, ierr

   IF (.NOT.ALLOCATED(ibufer)) THEN
      ALLOCATE(ibufer(0:N_box_Vx_e), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : SNAP_ELECTRON_2D_VDF_WALLS : Error in ALLOCATE ibufer !!!")', Rank_of_process
         PRINT '(2x,"The program will be terminated now :(")'
         STOP
      END IF
   END IF

   IF (.NOT.ALLOCATED(ibufer2)) THEN
      ALLOCATE(ibufer2(0:N_box_Vx_e), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : SNAP_ELECTRON_2D_VDF_WALLS : Error in ALLOCATE ibufer2 !!!")', Rank_of_process
         PRINT '(2x,"The program will be terminated now :(")'
         STOP
      END IF
   END IF

   IF (Rank_of_process.NE.0) THEN   ! client process
! transmit data to the server
! ep_2vdf_rw
      DO j = 0, N_box_Vyz_e
         ibufer(0:N_box_Vx_e) = ep_2vdf_rw(0:N_box_Vx_e, j)
         ibufer2 = 0
         CALL MPI_REDUCE(ibufer, ibufer2, N_box_Vx_e + 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO
! ep_2vdf_lw
      DO j = 0, N_box_Vyz_e
         ibufer(0:N_box_Vx_e) = ep_2vdf_lw(0:N_box_Vx_e, j)
         ibufer2 = 0
         CALL MPI_REDUCE(ibufer, ibufer2, N_box_Vx_e + 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO
! es_2vdf_rw
      DO j = 0, N_box_Vyz_e
         ibufer(0:N_box_Vx_e) = es_2vdf_rw(0:N_box_Vx_e, j)
         ibufer2 = 0
         CALL MPI_REDUCE(ibufer, ibufer2, N_box_Vx_e + 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO
! es_2vdf_lw
      DO j = 0, N_box_Vyz_e
         ibufer(0:N_box_Vx_e) = es_2vdf_lw(0:N_box_Vx_e, j)
         ibufer2 = 0
         CALL MPI_REDUCE(ibufer, ibufer2, N_box_Vx_e + 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO

   ELSE                             ! server process
! receive data from clients
! ep_2vdf_rw
      DO j = 0, N_box_Vyz_e
         ibufer = 0
         ibufer2 = 0
         CALL MPI_REDUCE(ibufer2, ibufer, N_box_Vx_e + 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         ep_2vdf_rw(0:N_box_Vx_e, j) = ibufer(0:N_box_Vx_e)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO
! ep_2vdf_lw
      DO j = 0, N_box_Vyz_e
         ibufer = 0
         ibufer2 = 0
         CALL MPI_REDUCE(ibufer2, ibufer, N_box_Vx_e + 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         ep_2vdf_lw(0:N_box_Vx_e, j) = ibufer(0:N_box_Vx_e)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO
! es_2vdf_rw
      DO j = 0, N_box_Vyz_e
         ibufer = 0
         ibufer2 = 0
         CALL MPI_REDUCE(ibufer2, ibufer, N_box_Vx_e + 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         es_2vdf_rw(0:N_box_Vx_e, j) = ibufer(0:N_box_Vx_e)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO
! es_2vdf_lw
      DO j = 0, N_box_Vyz_e
         ibufer = 0
         ibufer2 = 0
         CALL MPI_REDUCE(ibufer2, ibufer, N_box_Vx_e + 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         es_2vdf_lw(0:N_box_Vx_e, j) = ibufer(0:N_box_Vx_e)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO

      e2vdfw_filename = '_TTTT_e2vdfw.dat'
      e2vdfw_filename(2:5) = snapnumber_txt

!!     OPEN (99, FILE = e2vdfw_filename)
!!     WRITE (99, '("# col 1 is |Vx| [in units of V_th_e = ",e12.5," m/s]")') V_Te_ms
!!     WRITE (99, '("# col 2 is |V_perp|=sqrt(Vy**2+Vz**2) [V_th_e]")')
!!     WRITE (99, '("#----------")')
!!     WRITE (99, '("# col 3 is VDF of electrons hitting      the left wall")')
!!     WRITE (99, '("# col 4 is VDF of electrons emitted from the left wall")')
!!     WRITE (99, '("# col 5 is VDF of electrons hitting      the right wall")')
!!     WRITE (99, '("# col 6 is VDF of electrons emitted from the right wall")')
!!     WRITE (99, '("#----------")')
!!     DO j = 0, N_box_Vyz_e
!!        DO i = 0, N_box_Vx_e
!!           WRITE (99, '(2(1x,f10.5),4(1x,i7))') &
!!                & evx_mid_of_box(i+1), &
!!                & evyz_mid_of_box(j+1), &
!!                & ep_2vdf_lw(i, j), &
!!                & es_2vdf_lw(i, j), &
!!                & ep_2vdf_rw(i, j), &
!!                & es_2vdf_rw(i, j)
!!        END DO
!!        WRITE  (99, '(" ")')
!!     END DO
!!     CLOSE (99, STATUS = 'KEEP')

   END IF

   ep_2vdf_lw = 0
   ep_2vdf_rw = 0
   es_2vdf_lw = 0
   es_2vdf_rw = 0

   IF (ALLOCATED(ibufer)) THEN
      DEALLOCATE(ibufer, STAT = DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : SNAP_ELECTRON_2D_VDF_WALLS : Error in DEALLOCATE ibufer !!!")', Rank_of_process
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(ibufer2)) THEN
      DEALLOCATE(ibufer2, STAT = DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : SNAP_ELECTRON_2D_VDF_WALLS : Error in DEALLOCATE ibufer2 !!!")', Rank_of_process
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

END SUBROUTINE SNAP_ELECTRON_2D_VDF_WALLS

!------------------------------------------------------
! produces the local ion velocity distribution functions
SUBROUTINE SNAP_LOCAL_ELECTRON_VDFS

   use mpi
   USE ParallelOperationValues
   USE Snapshots
   USE CurrentProblemValues

   IMPLICIT NONE

!  INCLUDE 'mpif.h'

   INTEGER k               ! particle index
   INTEGER indx_x          ! index of x-velocity box
   INTEGER indx_y          ! index of y-velocity box
   INTEGER indx_z          ! index of z-velocity box
   INTEGER indx_par        ! index of parallel-velocity box
   INTEGER loc             ! index of current location
   INTEGER cell            ! index of cell where the particle is

   CHARACTER(15) evxdf_filename
   CHARACTER(15) evydf_filename
   CHARACTER(15) evzdf_filename

   INTEGER N_in_loc(1:N_of_all_vdf_locs)     ! number of particles inside the regions, where the distribution function is calculated
   ! these values are used for normalization

   REAL(8) temp_arr(1:(3*N_of_all_vdf_locs))

   INTEGER, ALLOCATABLE :: ibufer(:)
   INTEGER, ALLOCATABLE :: ibufer2(:)
   INTEGER ALLOC_ERR, DEALLOC_ERR, ierr

!  N_box_Vx_e      = Ve_x_max * N_box_vel - 1
! -N_box_Vx_e      = N_box_Vx_e_low  this will save three operations of addition and three operations of sign changing for each particle
!  N_box_Vx_e + 1  = N_box_Vx_e_top
! -N_box_Vyz_e     = N_box_Vyz_e_low
!  N_box_Vyz_e + 1 = N_box_Vyz_e_top

   IF (N_of_all_vdf_locs.LT.1) RETURN     ! quit if creation of local distributions was not requested

   ALLOCATE(ibufer(1:(2*MAX(N_box_Vx_e_top,N_box_Vyz_e_top))), STAT=ALLOC_ERR)   ! assume that 2*MAX(N_box_Vx_e_top,N_box_Vyz_e_top) > N_of_all_vdf_locs
   IF(ALLOC_ERR.NE.0)THEN
      PRINT '(2x,"Process ",i3," : SNAP_LOCAL_ELECTRON_VDFS : Error in ALLOCATE ibufer !!!")', Rank_of_process
      PRINT '(2x,"The program will be terminated now :(")'
      STOP
   END IF

   ALLOCATE(ibufer2(1:(2*MAX(N_box_Vx_e_top,N_box_Vyz_e_top))), STAT=ALLOC_ERR)   ! assume that 2*MAX(N_box_Vx_e_top,N_box_Vyz_e_top) > N_of_all_vdf_locs
   IF(ALLOC_ERR.NE.0)THEN
      PRINT '(2x,"Process ",i3," : SNAP_LOCAL_ELECTRON_VDFS : Error in ALLOCATE ibufer2 !!!")', Rank_of_process
      PRINT '(2x,"The program will be terminated now :(")'
      STOP
   END IF

   IF (Rank_of_process.GT.0) THEN                                                                             ! client >>>

      N_in_loc    = 0

! calculate arrays of distribution function
      DO k = 1, N_part(1)                             !

         if (flag_evxdf) then
            indx_x = INT(VX_of_spec(1)%part(k))
            IF (VX_of_spec(1)%part(k).GT.0.0_8) indx_x = indx_x + 1
            IF ((indx_x.LT.N_box_Vx_e_low).OR.(indx_x.GT.N_box_Vx_e_top))   CYCLE ! skip electron if it is too fast
         end if

         if (flag_evydf) then
            indx_y = INT(VY_of_spec(1)%part(k))
            IF (VY_of_spec(1)%part(k).GT.0.0_8) indx_y = indx_y + 1
            IF ((indx_y.LT.N_box_Vyz_e_low).OR.(indx_y.GT.N_box_Vyz_e_top)) CYCLE ! skip electron if it is too fast
         end if

         if (flag_evzdf) then
            indx_z = INT(VZ_of_spec(1)%part(k))
            IF (VZ_of_spec(1)%part(k).GT.0.0_8) indx_z = indx_z + 1
            IF ((indx_z.LT.N_box_Vyz_e_low).OR.(indx_z.GT.N_box_Vyz_e_top)) CYCLE ! skip electron if it is too fast
         end if

         cell = INT(X_of_spec(1)%part(k))

         DO loc = 1, N_of_all_vdf_locs
            IF (cell.LE.Vdf_location_bnd(loc)) THEN

               if (flag_evxdf) e_vxdf_loc(indx_x, loc) = e_vxdf_loc(indx_x, loc) + 1
               if (flag_evydf) e_vydf_loc(indx_y, loc) = e_vydf_loc(indx_y, loc) + 1
               if (flag_evzdf) e_vzdf_loc(indx_z, loc) = e_vzdf_loc(indx_z, loc) + 1
               N_in_loc(loc) = N_in_loc(loc) + 1

               IF (Tag_of_spec(1)%part(k).EQ.eTag_Emit_Left) THEN
                  ebl_vxdf_loc(indx_x, loc) = ebl_vxdf_loc(indx_x, loc) + 1
                  ebl_vydf_loc(indx_y, loc) = ebl_vydf_loc(indx_y, loc) + 1
                  ebl_vzdf_loc(indx_z, loc) = ebl_vzdf_loc(indx_z, loc) + 1
               ELSE IF (Tag_of_spec(1)%part(k).EQ.eTag_Emit_Right) THEN
                  ebr_vxdf_loc(indx_x, loc) = ebr_vxdf_loc(indx_x, loc) + 1
                  ebr_vydf_loc(indx_y, loc) = ebr_vydf_loc(indx_y, loc) + 1
                  ebr_vzdf_loc(indx_z, loc) = ebr_vzdf_loc(indx_z, loc) + 1
               END IF

               EXIT
            END IF
         END DO

      END DO

! transmit data to the server
      ! N_in_loc
      ibufer(1:N_of_all_vdf_locs)  = N_in_loc(1:N_of_all_vdf_locs)
      ibufer2(1:N_of_all_vdf_locs) = 0
      CALL MPI_REDUCE(ibufer, ibufer2, N_of_all_vdf_locs, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

      ! distributions . . .
      DO loc = 1, N_of_all_vdf_locs          ! cycle over locations
      ! common
         if (flag_evxdf) then
            ! e_vxdf_loc
            ibufer(1:(N_box_Vx_e_top+N_box_Vx_e_top))  = e_vxdf_loc(N_box_Vx_e_low:N_box_Vx_e_top, loc)
            ibufer2(1:(N_box_Vx_e_top+N_box_Vx_e_top)) = 0
            CALL MPI_REDUCE(ibufer, ibufer2, (N_box_Vx_e_top+N_box_Vx_e_top), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
         end if

         if (flag_evydf) then
            ! e_vydf_loc
            ibufer(1:(N_box_Vyz_e_top+N_box_Vyz_e_top))  = e_vydf_loc(N_box_Vyz_e_low:N_box_Vyz_e_top, loc)
            ibufer2(1:(N_box_Vyz_e_top+N_box_Vyz_e_top)) = 0
            CALL MPI_REDUCE(ibufer, ibufer2, (N_box_Vyz_e_top+N_box_Vyz_e_top), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
         end if

         if (flag_evzdf) then
            ! e_vzdf_loc
            ibufer(1:(N_box_Vyz_e_top+N_box_Vyz_e_top))  = e_vzdf_loc(N_box_Vyz_e_low:N_box_Vyz_e_top, loc)
            ibufer2(1:(N_box_Vyz_e_top+N_box_Vyz_e_top)) = 0
            CALL MPI_REDUCE(ibufer, ibufer2, (N_box_Vyz_e_top+N_box_Vyz_e_top), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
         end if

      ! emitted from the left
! ebl_vxdf_loc
         ibufer(1:(N_box_Vx_e_top+N_box_Vx_e_top))  = ebl_vxdf_loc(N_box_Vx_e_low:N_box_Vx_e_top, loc)
         ibufer2(1:(N_box_Vx_e_top+N_box_Vx_e_top)) = 0
         CALL MPI_REDUCE(ibufer, ibufer2, (N_box_Vx_e_top+N_box_Vx_e_top), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
! ebl_vydf_loc
         ibufer(1:(N_box_Vyz_e_top+N_box_Vyz_e_top))  = ebl_vydf_loc(N_box_Vyz_e_low:N_box_Vyz_e_top, loc)
         ibufer2(1:(N_box_Vyz_e_top+N_box_Vyz_e_top)) = 0
         CALL MPI_REDUCE(ibufer, ibufer2, (N_box_Vyz_e_top+N_box_Vyz_e_top), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
! ebl_vzdf_loc
         ibufer(1:(N_box_Vyz_e_top+N_box_Vyz_e_top))  = ebl_vzdf_loc(N_box_Vyz_e_low:N_box_Vyz_e_top, loc)
         ibufer2(1:(N_box_Vyz_e_top+N_box_Vyz_e_top)) = 0
         CALL MPI_REDUCE(ibufer, ibufer2, (N_box_Vyz_e_top+N_box_Vyz_e_top), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

      ! emitted from the right
! ebr_vxdf_loc
         ibufer(1:(N_box_Vx_e_top+N_box_Vx_e_top))  = ebr_vxdf_loc(N_box_Vx_e_low:N_box_Vx_e_top, loc)
         ibufer2(1:(N_box_Vx_e_top+N_box_Vx_e_top)) = 0
         CALL MPI_REDUCE(ibufer, ibufer2, (N_box_Vx_e_top+N_box_Vx_e_top), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
! ebr_vydf_loc
         ibufer(1:(N_box_Vyz_e_top+N_box_Vyz_e_top))  = ebr_vydf_loc(N_box_Vyz_e_low:N_box_Vyz_e_top, loc)
         ibufer2(1:(N_box_Vyz_e_top+N_box_Vyz_e_top)) = 0
         CALL MPI_REDUCE(ibufer, ibufer2, (N_box_Vyz_e_top+N_box_Vyz_e_top), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
! ebr_vzdf_loc
         ibufer(1:(N_box_Vyz_e_top+N_box_Vyz_e_top))  = ebr_vzdf_loc(N_box_Vyz_e_low:N_box_Vyz_e_top, loc)
         ibufer2(1:(N_box_Vyz_e_top+N_box_Vyz_e_top)) = 0
         CALL MPI_REDUCE(ibufer, ibufer2, (N_box_Vyz_e_top+N_box_Vyz_e_top), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO                                 ! end of cycle over locations

   ELSE                                                                                     ! server >>>

! receive data from clients
! N_in_loc
      ibufer  = 0
      ibufer2 = 0
      CALL MPI_REDUCE(ibufer2, ibufer, N_of_all_vdf_locs, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      N_in_loc(1:N_of_all_vdf_locs) = ibufer(1:N_of_all_vdf_locs)

      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

! distributions . . .
      DO loc = 1, N_of_all_vdf_locs          ! cycle over locations
      ! common
         if (flag_evxdf) then
            ! e_vxdf_loc
            ibufer  = 0
            ibufer2 = 0
            CALL MPI_REDUCE(ibufer2, ibufer, (N_box_Vx_e_top+N_box_Vx_e_top), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            e_vxdf_loc(N_box_Vx_e_low:N_box_Vx_e_top, loc) = ibufer(1:(N_box_Vx_e_top+N_box_Vx_e_top))
            CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
         end if

         if (flag_evydf) then
            ! e_vydf_loc
            ibufer = 0
            ibufer2 = 0
            CALL MPI_REDUCE(ibufer2, ibufer, (N_box_Vyz_e_top+N_box_Vyz_e_top), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            e_vydf_loc(N_box_Vyz_e_low:N_box_Vyz_e_top, loc) = ibufer(1:(N_box_Vyz_e_top+N_box_Vyz_e_top))
            CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
         end if

         if (flag_evzdf) then
            ! e_vzdf_loc
            ibufer = 0
            ibufer2 = 0
            CALL MPI_REDUCE(ibufer2, ibufer, (N_box_Vyz_e_top+N_box_Vyz_e_top), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            e_vzdf_loc(N_box_Vyz_e_low:N_box_Vyz_e_top, loc) = ibufer(1:(N_box_Vyz_e_top+N_box_Vyz_e_top))
            CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
         end if

      ! emitted from the left
! ebl_vxdf_loc
         ibufer  = 0
         ibufer2 = 0
         CALL MPI_REDUCE(ibufer2, ibufer, (N_box_Vx_e_top+N_box_Vx_e_top), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         ebl_vxdf_loc(N_box_Vx_e_low:N_box_Vx_e_top, loc) = ibufer(1:(N_box_Vx_e_top+N_box_Vx_e_top))
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
! ebl_vydf_loc
         ibufer = 0
         ibufer2 = 0
         CALL MPI_REDUCE(ibufer2, ibufer, (N_box_Vyz_e_top+N_box_Vyz_e_top), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         ebl_vydf_loc(N_box_Vyz_e_low:N_box_Vyz_e_top, loc) = ibufer(1:(N_box_Vyz_e_top+N_box_Vyz_e_top))
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
! ebl_vzdf_loc
         ibufer = 0
         ibufer2 = 0
         CALL MPI_REDUCE(ibufer2, ibufer, (N_box_Vyz_e_top+N_box_Vyz_e_top), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         ebl_vzdf_loc(N_box_Vyz_e_low:N_box_Vyz_e_top, loc) = ibufer(1:(N_box_Vyz_e_top+N_box_Vyz_e_top))
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

      ! emitted from the right
! ebr_vxdf_loc
         ibufer  = 0
         ibufer2 = 0
         CALL MPI_REDUCE(ibufer2, ibufer, (N_box_Vx_e_top+N_box_Vx_e_top), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         ebr_vxdf_loc(N_box_Vx_e_low:N_box_Vx_e_top, loc) = ibufer(1:(N_box_Vx_e_top+N_box_Vx_e_top))
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
! ebr_vydf_loc
         ibufer = 0
         ibufer2 = 0
         CALL MPI_REDUCE(ibufer2, ibufer, (N_box_Vyz_e_top+N_box_Vyz_e_top), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         ebr_vydf_loc(N_box_Vyz_e_low:N_box_Vyz_e_top, loc) = ibufer(1:(N_box_Vyz_e_top+N_box_Vyz_e_top))
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
! ebr_vzdf_loc
         ibufer = 0
         ibufer2 = 0
         CALL MPI_REDUCE(ibufer2, ibufer, (N_box_Vyz_e_top+N_box_Vyz_e_top), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         ebr_vzdf_loc(N_box_Vyz_e_low:N_box_Vyz_e_top, loc) = ibufer(1:(N_box_Vyz_e_top+N_box_Vyz_e_top))
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO                                 ! end of cycle over locations

      if (flag_evxdf) then
         ! EVDF(VX)
         evxdf_filename = '_TTTT_evxdf.dat'
         evxdf_filename(2:5) = snapnumber_txt
         OPEN (99, FILE = evxdf_filename)
         WRITE (99, '("# col   1 is Vx [in units of V_th_e = ",e12.5," m/s]")') V_Te_ms
         WRITE (99, '("#----------")')
         DO loc = 1, N_of_all_vdf_locs
            WRITE (99, '("# location ",i2," with ",i8," macroparticles")') loc, N_in_loc(loc)
            WRITE (99, '("# col ",i3," is VDF(Vx) of ALL electrons")')                         3*loc-1 !1+1+3*(k-1)
            WRITE (99, '("# col ",i3," is VDF(Vx) of electrons emitted from the LEFT wall")')  3*loc   !1+2+3*(k-1)
            WRITE (99, '("# col ",i3," is VDF(Vx) of electrons emitted from the RIGHT wall")') 3*loc+1 !1+3+3*(k-1)
            WRITE (99, '("#----------")')
         END DO
         DO indx_x = N_box_Vx_e_low, N_box_Vx_e_top
            DO loc = 1, N_of_all_vdf_locs
               temp_arr(3*loc-2) = DBLE(  e_vxdf_loc(indx_x,loc)) * DBLE(N_box_vel)
               temp_arr(3*loc-1) = DBLE(ebl_vxdf_loc(indx_x,loc)) * DBLE(N_box_vel)
               temp_arr(3*loc)   = DBLE(ebr_vxdf_loc(indx_x,loc)) * DBLE(N_box_vel)
            END DO
            WRITE (99, '(2x,f10.5,297(1x,e12.5))') evx_mid_of_box(indx_x), temp_arr(1:(3*N_of_all_vdf_locs))
         END DO
         CLOSE (99, STATUS = 'KEEP')
      end if

      if (flag_evydf) then
         ! EVDF(VY)
         evydf_filename = '_TTTT_evydf.dat'
         evydf_filename(2:5) = snapnumber_txt
         OPEN (99, FILE = evydf_filename)
         WRITE (99, '("# col   1 is Vy [in units of V_th_e = ",e12.5," m/s]")') V_Te_ms
         WRITE (99, '("#----------")')
         DO loc = 1, N_of_all_vdf_locs
            WRITE (99, '("# location ",i2," with ",i8," macroparticles")') loc, N_in_loc(loc)
            WRITE (99, '("# col ",i3," is VDF(Vy) of ALL electrons")')                         3*loc-1 !1+1+3*(k-1)
            WRITE (99, '("# col ",i3," is VDF(Vy) of electrons emitted from the LEFT wall")')  3*loc   !1+2+3*(k-1)
            WRITE (99, '("# col ",i3," is VDF(Vy) of electrons emitted from the RIGHT wall")') 3*loc+1 !1+3+3*(k-1)
            WRITE (99, '("#----------")')
         END DO
         DO indx_y = N_box_Vyz_e_low, N_box_Vyz_e_top
            DO loc = 1, N_of_all_vdf_locs
               temp_arr(3*loc-2) = DBLE(  e_vydf_loc(indx_y,loc)) * DBLE(N_box_vel)
               temp_arr(3*loc-1) = DBLE(ebl_vydf_loc(indx_y,loc)) * DBLE(N_box_vel)
               temp_arr(3*loc)   = DBLE(ebr_vydf_loc(indx_y,loc)) * DBLE(N_box_vel)
            END DO
            WRITE (99, '(2x,f10.5,297(1x,e12.5))') evyz_mid_of_box(indx_y), temp_arr(1:(3*N_of_all_vdf_locs))
         END DO
         CLOSE (99, STATUS = 'KEEP')
      end if

      if (flag_evzdf) then
         ! EVDF(VZ)
         evzdf_filename = '_TTTT_evzdf.dat'
         evzdf_filename(2:5) = snapnumber_txt
         OPEN (99, FILE = evzdf_filename)
         WRITE (99, '("# col   1 is Vz [in units of V_th_e = ",e12.5," m/s]")') V_Te_ms
         WRITE (99, '("#----------")')
         DO loc = 1, N_of_all_vdf_locs
            WRITE (99, '("# location ",i2," with ",i8," macroparticles")') loc, N_in_loc(loc)
            WRITE (99, '("# col ",i3," is VDF(Vz) of ALL electrons")')                         3*loc-1 !1+1+3*(k-1)
            WRITE (99, '("# col ",i3," is VDF(Vz) of electrons emitted from the LEFT wall")')  3*loc   !1+2+3*(k-1)
            WRITE (99, '("# col ",i3," is VDF(Vz) of electrons emitted from the RIGHT wall")') 3*loc+1 !1+3+3*(k-1)
            WRITE (99, '("#----------")')
         END DO
         DO indx_z = N_box_Vyz_e_low, N_box_Vyz_e_top
            DO loc = 1, N_of_all_vdf_locs
               temp_arr(3*loc-2)= DBLE(  e_vzdf_loc(indx_z,loc)) * DBLE(N_box_vel)
               temp_arr(3*loc-1)= DBLE(ebl_vzdf_loc(indx_z,loc)) * DBLE(N_box_vel)
               temp_arr(3*loc)  = DBLE(ebr_vzdf_loc(indx_z,loc)) * DBLE(N_box_vel)
            END DO
            WRITE (99, '(2x,f10.5,297(1x,e12.5))') evyz_mid_of_box(indx_z), temp_arr(1:(3*N_of_all_vdf_locs))
         END DO
         CLOSE (99, STATUS = 'KEEP')
      end if

   END IF

   
   if (flag_evxdf) e_vxdf_loc = 0
   if (flag_evydf) e_vydf_loc = 0
   if (flag_evzdf) e_vzdf_loc = 0

   ebl_vxdf_loc = 0
   ebl_vydf_loc = 0
   ebl_vzdf_loc = 0

   ebr_vxdf_loc = 0
   ebr_vydf_loc = 0
   ebr_vzdf_loc = 0

   IF (ALLOCATED(ibufer)) THEN
      DEALLOCATE(ibufer, STAT = DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : SNAP_LOCAL_ELECTRON_VDFS : Error in DEALLOCATE ibufer !!!")', Rank_of_process
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(ibufer2)) THEN
      DEALLOCATE(ibufer2, STAT = DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : SNAP_LOCAL_ELECTRON_VDFS : Error in DEALLOCATE ibufer2 !!!")', Rank_of_process
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

END SUBROUTINE SNAP_LOCAL_ELECTRON_VDFS

!! EVDF(VX,VY) and EVDF(VX,VZ) --------------------------- omit this part for now due to the large file size
!---------------   it is easier to save the complete phase plane and then restore the 2D EVDFs at any location one wants
!
!        e2vxy2vxzdf_filename = '_TTTT_e2vxy2vxzdf.dat'
!        e2vxy2vxzdf_filename(2:5) = snapnumber_txt
!        OPEN (99, FILE = e2vxvydf_filename)
!        WRITE (99, '("# col   1 is Vx [in units of V_th_e = ",e12.5," m/s]")') V_Te_ms
!        WRITE (99, '("# col   2 is Vy or Vz [V_th_e]")')
!        WRITE (99, '("#----------")')
!        DO k = 1, N_of_all_vdf_locs
!           WRITE (99, '("# location ",i2)') k
!           WRITE (99, '("# col ",i3," is VDF(Vx,Vy) of all electrons")')                         6*k-3 !2+1+6*(k-1)
!           WRITE (99, '("# col ",i3," is VDF(Vx,Vy) of electrons emitted from the left wall")')  6*k-2 !2+2+6*(k-1)
!           WRITE (99, '("# col ",i3," is VDF(Vx,Vy) of electrons emitted from the right wall")') 6*k-1 !2+3+6*(k-1)
!           WRITE (99, '("# col ",i3," is VDF(Vx,Vz) of all electrons")')                         6*k   !2+4+6*(k-1)
!           WRITE (99, '("# col ",i3," is VDF(Vx,Vz) of electrons emitted from the left wall")')  6*k+1 !2+5+6*(k-1)
!           WRITE (99, '("# col ",i3," is VDF(Vx,Vz) of electrons emitted from the right wall")') 6*k+2 !2+6+6*(k-1)
!           WRITE (99, '("#----------")')
!        END DO
!        DO j = N_box_Vyz_e_low, N_box_Vyz_e_top
!           DO i = N_box_Vx_e_low, N_box_Vx_e_top
!              DO k = 1, 6*N_of_all_vdf_locs-5, 6
!                 temp_arr(k)   =   e_2vxvydf_loc(i,j,k)
!                 temp_arr(k+1) = ebl_2vxvydf_loc(i,j,k)
!                 temp_arr(k+2) = ebr_2vxvydf_loc(i,j,k)
!                 temp_arr(k+3) =   e_2vxvzdf_loc(i,j,k)
!                 temp_arr(k+4) = ebl_2vxvzdf_loc(i,j,k)
!                 temp_arr(k+5) = ebr_2vxvzdf_loc(i,j,k)
!              END DO
!              WRITE (99, '(2(1x,f10.5),594(1x,i7))') &    ! 99*6=594
!                   & evx_mid_of_box(i), &
!                   & evyz_mid_of_box(j), &
!                   & temp_arr(1:(6*N_of_all_vdf_locs))
!           END DO
!           WRITE (99, '(" ")')
!        END DO
!        CLOSE (99, STATUS = 'KEEP')

!------------------------------------------
! produces the phase planes for electrons
SUBROUTINE SNAP_ELECTRON_PHASE_PLANES

   use mpi
   USE ParallelOperationValues
   USE Snapshots
   USE CurrentProblemValues

   IMPLICIT NONE

!  INCLUDE 'mpif.h'

   REAL(8), ALLOCATABLE :: dbufer(:)
   INTEGER N_of_beam_left
   INTEGER N_of_beam_right
   INTEGER bufer_len, transm_len
   INTEGER ALLOC_ERR, DEALLOC_ERR, ierr
   INTEGER stattus(MPI_STATUS_SIZE)

   INTEGER k               ! particle index
   INTEGER n               ! index in the bufer of transmission
   INTEGER skip_count      ! counts particles to skip in case of beam/emitted particles

   CHARACTER(13) epp_filename
   CHARACTER(15) elwpp_filename
   CHARACTER(15) erwpp_filename

   IF (Flag_pp_creation.EQ.0) RETURN

!     PRINT '(2x,"Process ",i3," : Creating electron phase planes...")', Rank_of_process

   IF (Rank_of_process.GT.0) THEN                                                                             ! client >>>

! determine the number of beam particles in the process
      N_of_beam_left = 0
      N_of_beam_right = 0
      DO k = 1, N_part(1)
         IF (Tag_of_spec(1)%part(k).EQ.eTag_Emit_Left)  N_of_beam_left  = N_of_beam_left  + 1
         IF (Tag_of_spec(1)%part(k).EQ.eTag_Emit_Right) N_of_beam_right = N_of_beam_right + 1
      END DO

! define the maximal necessary buffer length
!     bufer_len = 4 * MAX(N_of_beam_left, N_of_beam_right, ((N_part(1) - 1) / (N_to_skip + 1) + 1))

      bufer_len = 4 * MAX( (N_of_beam_left-1)  / (N_to_skip_left  + 1) + 10, &
      & (N_of_beam_right-1) / (N_to_skip_right + 1) + 10, &
      & (N_part(1)-1)       / (N_to_skip       + 1) + 10  )

! send the buffer length to the server process
      CALL MPI_SEND(bufer_len, 1, MPI_INTEGER, 0, 10, MPI_COMM_WORLD, ierr)

! allocate the buffer array
      IF (.NOT.ALLOCATED(dbufer)) THEN
         ALLOCATE(dbufer(1:bufer_len), STAT=ALLOC_ERR)
         IF(ALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : SNAP_ELECTRON_PHASE_PLANES : Error in ALLOCATE dbufer !!!")', Rank_of_process
            PRINT '(2x,"The program will be terminated now :(")'
            STOP
         END IF
      END IF

! for the non-beam particles only (but note that the collided particles are included here)
      n = 1
      DO k = 1, N_part(1), N_to_skip + 1                          ! fill the buffer with the particles data
         IF ((Tag_of_spec(1)%part(k).NE.eTag_Emit_Left).AND.(Tag_of_spec(1)%part(k).NE.eTag_Emit_Right)) THEN
            dbufer(n)   =  X_of_spec(1)%part(k)
            dbufer(n+1) = VX_of_spec(1)%part(k)
            dbufer(n+2) = VY_of_spec(1)%part(k)
            dbufer(n+3) = VZ_of_spec(1)%part(k)
            n = n + 4
         END IF
      END DO
      transm_len = n - 1      ! now transm_len is the length of the transmission

      IF (transm_len.GT.bufer_len) THEN
         PRINT '(2x,"Process ",i3," : Error in SNAP_ELECTRON_PHASE_PLANES !!!")', Rank_of_process
         PRINT '(2x,"Index in the bufer ",i8," exceeds the threshold value ",i8)', transm_len, bufer_len
         PRINT '(2x,"The program will be terminated now :(")'
         STOP
      END IF

! transmit data
! bufer_len
      CALL MPI_SEND(transm_len, 1, MPI_INTEGER, 0, 11, MPI_COMM_WORLD, ierr)
! dbufer, if the length is non-zero
      IF (transm_len.GT.0) THEN
         CALL MPI_SEND(dbufer, transm_len, MPI_DOUBLE_PRECISION, 0, 12, MPI_COMM_WORLD, ierr)
      END IF

! for electrons emitted from the left wall
      n = 1
      skip_count = N_to_skip_left !0
      DO k = 1, N_part(1)                                ! fill the buffer with the particles data
         IF (Tag_of_spec(1)%part(k).EQ.eTag_Emit_Left) THEN
            skip_count = skip_count+1
            IF (skip_count.GT.N_to_skip_left) THEN
               dbufer(n)   =  X_of_spec(1)%part(k)
               dbufer(n+1) = VX_of_spec(1)%part(k)
               dbufer(n+2) = VY_of_spec(1)%part(k)
               dbufer(n+3) = VZ_of_spec(1)%part(k)
               n = n + 4
               skip_count = 0
            END IF
         END IF
      END DO
      transm_len = n - 1      ! now transm_len is the length of transmission

      IF (transm_len.GT.bufer_len) THEN
         PRINT '(2x,"Process ",i3," : Error in SNAP_ELECTRON_PHASE_PLANES !!!")', Rank_of_process
         PRINT '(2x,"last index in dbufer ",i8," exceeds the limit for particles emitted from the left wall ",i8)', transm_len, bufer_len
         PRINT '(2x,"The program will be terminated now :(")'
         STOP
      END IF

! transmit data
! bufer_len
      CALL MPI_SEND(transm_len, 1, MPI_INTEGER, 0, 13, MPI_COMM_WORLD, ierr)
! dbufer, if the length is non-zero
      IF (transm_len.GT.0) THEN
         CALL MPI_SEND(dbufer, transm_len, MPI_DOUBLE_PRECISION, 0, 14, MPI_COMM_WORLD, ierr)
      END IF

! for electrons emitted from the right wall
      n = 1
      skip_count = N_to_skip_right !0
      DO k = 1, N_part(1)                                ! fill the buffer with the particles data
         IF (Tag_of_spec(1)%part(k).EQ.eTag_Emit_Right) THEN
            skip_count = skip_count+1
            IF (skip_count.GT.N_to_skip_right) THEN
               dbufer(n)   =  X_of_spec(1)%part(k)
               dbufer(n+1) = VX_of_spec(1)%part(k)
               dbufer(n+2) = VY_of_spec(1)%part(k)
               dbufer(n+3) = VZ_of_spec(1)%part(k)
               n = n + 4
               skip_count = 0
            END IF
         END IF
      END DO
      transm_len = n - 1

      IF (transm_len.GT.bufer_len) THEN
         PRINT '(2x,"Process ",i3," : Error in SNAP_ELECTRON_PHASE_PLANES !!!")', Rank_of_process
         PRINT '(2x,"last index in dbufer ",i8," exceeds the limit for particles emitted from the right wall ",i8)', transm_len, bufer_len
         PRINT '(2x,"The program will be terminated now :(")'
         STOP
      END IF

! transmit data
! bufer_len
      CALL MPI_SEND(transm_len, 1, MPI_INTEGER, 0, 15, MPI_COMM_WORLD, ierr)
! dbufer, if the length is non-zero
      IF (transm_len.GT.0) THEN
         CALL MPI_SEND(dbufer, transm_len, MPI_DOUBLE_PRECISION, 0, 16, MPI_COMM_WORLD, ierr)
      END IF

! deallocate the buffer
      IF (ALLOCATED(dbufer)) THEN
         DEALLOCATE(dbufer, STAT = DEALLOC_ERR)
         IF(DEALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : SNAP_ELECTRON_PHASE_PLANES : Error in DEALLOCATE dbufer !!!")', Rank_of_process
            PRINT *, 'The program will be terminated now :('
            STOP
         END IF
      END IF

   ELSE                                                                           ! server >>>>

      PRINT '(2x,"Process ",i3," : Creating electron phase planes...")', Rank_of_process

! Produce new filenames
      epp_filename      =  '_TTTT_epp.dat'
      epp_filename(2:5) = snapnumber_txt
      elwpp_filename      =  '_TTTT_elwpp.dat'
      elwpp_filename(2:5) = snapnumber_txt
      erwpp_filename      =  '_TTTT_erwpp.dat'
      erwpp_filename(2:5) = snapnumber_txt
      OPEN (41, FILE =   epp_filename)
      OPEN (42, FILE = elwpp_filename)
      OPEN (43, FILE = erwpp_filename)

      DO k = 1, N_of_processes - 1

! get the maximal buffer length
         CALL MPI_PROBE(k, 10, MPI_COMM_WORLD, stattus, ierr)
         CALL MPI_RECV(bufer_len, 1, MPI_INTEGER, k, 10, MPI_COMM_WORLD, stattus, ierr)

! allocate the buffer
         IF (.NOT.ALLOCATED(dbufer)) THEN
            ALLOCATE(dbufer(1:bufer_len), STAT=ALLOC_ERR)
            IF(ALLOC_ERR.NE.0)THEN
               PRINT '(2x,"Process ",i3," : SNAP_ELECTRON_PHASE_PLANES : Error in ALLOCATE dbufer !!!")', Rank_of_process
               PRINT '(2x,"The program will be terminated now :(")'
               STOP
            END IF
         END IF

! receive the non beam particles
         transm_len = 0
         CALL MPI_PROBE(k, 11, MPI_COMM_WORLD, stattus, ierr)
         CALL MPI_RECV(transm_len, 1, MPI_INTEGER, k, 11, MPI_COMM_WORLD, stattus, ierr)

         IF (transm_len.GT.0) THEN
            dbufer = 0.0_8
            CALL MPI_PROBE(k, 12, MPI_COMM_WORLD, stattus, ierr)
            CALL MPI_RECV(dbufer, transm_len, MPI_DOUBLE_PRECISION, k, 12, MPI_COMM_WORLD, stattus, ierr)
            DO n = 1, transm_len - 3, 4
               WRITE (41, '(4(2x,e12.5))') dbufer(n), dbufer(n+1), dbufer(n+2), dbufer(n+3)
            END DO
         END IF

! receive particles of the beam emitted from the left wall
         transm_len = 0
         CALL MPI_PROBE(k, 13, MPI_COMM_WORLD, stattus, ierr)
         CALL MPI_RECV(transm_len, 1, MPI_INTEGER, k, 13, MPI_COMM_WORLD, stattus, ierr)

         IF (transm_len.GT.0) THEN
            dbufer = 0.0_8
            CALL MPI_PROBE(k, 14, MPI_COMM_WORLD, stattus, ierr)
            CALL MPI_RECV(dbufer, transm_len, MPI_DOUBLE_PRECISION, k, 14, MPI_COMM_WORLD, stattus, ierr)
            DO n = 1, transm_len - 3, 4
               WRITE (42, '(4(2x,e12.5))') dbufer(n), dbufer(n+1), dbufer(n+2), dbufer(n+3)
            END DO
         END IF

! receive particles of the beam emitted from the right wall
         transm_len = 0
         CALL MPI_PROBE(k, 15, MPI_COMM_WORLD, stattus, ierr)
         CALL MPI_RECV(transm_len, 1, MPI_INTEGER, k, 15, MPI_COMM_WORLD, stattus, ierr)

         IF (transm_len.GT.0) THEN
            dbufer = 0.0_8
            CALL MPI_PROBE(k, 16, MPI_COMM_WORLD, stattus, ierr)
            CALL MPI_RECV(dbufer, transm_len, MPI_DOUBLE_PRECISION, k, 16, MPI_COMM_WORLD, stattus, ierr)
            DO n = 1, transm_len - 3, 4
               WRITE (43, '(4(2x,e12.5))') dbufer(n), dbufer(n+1), dbufer(n+2), dbufer(n+3)
            END DO
         END IF

! deallocate the buffer
         IF (ALLOCATED(dbufer)) THEN
            DEALLOCATE(dbufer, STAT = DEALLOC_ERR)
            IF(DEALLOC_ERR.NE.0)THEN
               PRINT '(2x,"Process ",i3," : SNAP_ELECTRON_PHASE_PLANES : Error in DEALLOCATE dbufer !!!")', Rank_of_process
               PRINT *, 'The program will be terminated now :('
               STOP
            END IF
         END IF

      END DO

      CLOSE (41, STATUS = 'KEEP')
      CLOSE (42, STATUS = 'KEEP')
      CLOSE (43, STATUS = 'KEEP')

   END IF

END SUBROUTINE SNAP_ELECTRON_PHASE_PLANES

!------------------------------------------------------
! produces the local ion velocity distribution functions
SUBROUTINE SNAP_LOCAL_ION_VDFS

   use mpi
   USE ParallelOperationValues
   USE Snapshots
   USE CurrentProblemValues

   IMPLICIT NONE

!  INCLUDE 'mpif.h'

   INTEGER k               ! particle index
   INTEGER indx_x          ! index of velocity box
   ! INTEGER indx_y          ! index of y-velocity box
   ! INTEGER indx_z          ! index of z-velocity box
   INTEGER loc             ! index of current location
   INTEGER cell            ! index of cell where the particle is

   CHARACTER(15) ivxdf_filename
   ! CHARACTER(15) ivydf_filename
   ! CHARACTER(15) ivzdf_filename

   INTEGER N_in_loc(1:N_of_all_vdf_locs)     ! number of particles inside the regions, where the distribution function is calculated
   ! these values are used for normalization
   REAL(8) temp_arr(1:(3*N_of_all_vdf_locs))

   INTEGER, ALLOCATABLE :: ibufer(:)
   INTEGER, ALLOCATABLE :: ibufer2(:)
   INTEGER ALLOC_ERR, DEALLOC_ERR, ierr

! N_box_Vx_i_low = -N_box_Vx_i                   ! this will save one addition and one sign changing for each ion
! N_box_Vx_i_top =  N_box_Vx_i + 1               !

   IF (N_of_all_vdf_locs.LT.1) RETURN     ! quit if creation of local distributions was not requested

   IF (.NOT.ALLOCATED(ibufer)) THEN
      ALLOCATE(ibufer(1:(2*N_box_Vx_i_top)), STAT=ALLOC_ERR)   ! assume that 2*N_box_Vx_i_top > N_of_all_vdf_locs
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : SNAP_LOCAL_ION_VDFS : Error in ALLOCATE ibufer !!!")', Rank_of_process
         PRINT '(2x,"The program will be terminated now :(")'
         STOP
      END IF
   END IF

   IF (.NOT.ALLOCATED(ibufer2)) THEN
      ALLOCATE(ibufer2(1:(2*N_box_Vx_i_top)), STAT=ALLOC_ERR)   ! assume that 2*N_box_Vx_i_top > N_of_all_vdf_locs
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : SNAP_LOCAL_ION_VDFS : Error in ALLOCATE ibufer2 !!!")', Rank_of_process
         PRINT '(2x,"The program will be terminated now :(")'
         STOP
      END IF
   END IF

   IF (Rank_of_process.GT.0) THEN                                                                             ! client >>>

      N_in_loc = 0

! calculate arrays of distribution function
      DO k = 1, N_part(2)                             !

         indx_x = INT( VX_of_spec(2)%part(k) * SQRT(Ms(2)) )
         IF (VX_of_spec(2)%part(k).GT.0.0_8) indx_x = indx_x + 1
         IF ((indx_x.LT.N_box_Vx_i_low).OR.(indx_x.GT.N_box_Vx_i_top)) CYCLE ! skip ions which are too fast

         ! indx_y = INT( VY_of_spec(2)%part(k) * SQRT(Ms(2)) )
         ! IF (VY_of_spec(2)%part(k).GT.0.0_8) indx_y = indx_y + 1
         ! IF ((indx_y.LT.N_box_Vx_i_low).OR.(indx_y.GT.N_box_Vx_i_top)) CYCLE ! skip ions which are too fast

         ! indx_z = INT( VZ_of_spec(2)%part(k) * SQRT(Ms(2)) )
         ! IF (VZ_of_spec(2)%part(k).GT.0.0_8) indx_z = indx_z + 1
         ! IF ((indx_z.LT.N_box_Vx_i_low).OR.(indx_z.GT.N_box_Vx_i_top)) CYCLE ! skip ions which are too fast

         cell = INT(X_of_spec(2)%part(k))
         DO loc = 1, N_of_all_vdf_locs
            IF (cell.LE.Vdf_location_bnd(loc)) THEN

               i_vxdf_loc(indx_x, loc) = i_vxdf_loc(indx_x, loc) + 1
               ! i_vydf_loc(indx_y, loc) = i_vydf_loc(indx_y, loc) + 1
               ! i_vzdf_loc(indx_z, loc) = i_vzdf_loc(indx_z, loc) + 1
               N_in_loc(loc) = N_in_loc(loc) + 1

               EXIT
            END IF
         END DO

      END DO
! transmit data to the server
! N_in_loc
      ibufer(1:N_of_all_vdf_locs)  = N_in_loc(1:N_of_all_vdf_locs)
      ibufer2(1:N_of_all_vdf_locs) = 0
      CALL MPI_REDUCE(ibufer, ibufer2, N_of_all_vdf_locs, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

! distributions . . .
      DO loc = 1, N_of_all_vdf_locs          ! cycle over locations
! i_vxdf_loc
         ibufer(1:(N_box_Vx_i_top+N_box_Vx_i_top))  = i_vxdf_loc(N_box_Vx_i_low:N_box_Vx_i_top, loc)
         ibufer2(1:(N_box_Vx_i_top+N_box_Vx_i_top)) = 0
         CALL MPI_REDUCE(ibufer, ibufer2, (N_box_Vx_i_top+N_box_Vx_i_top), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
! ! i_vydf_loc
!          ibufer(1:(N_box_Vx_i_top+N_box_Vx_i_top))  = i_vydf_loc(N_box_Vx_i_low:N_box_Vx_i_top, loc)
!          ibufer2(1:(N_box_Vx_i_top+N_box_Vx_i_top)) = 0
!          CALL MPI_REDUCE(ibufer, ibufer2, (N_box_Vx_i_top+N_box_Vx_i_top), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
!          CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
! ! i_vzdf_loc
!          ibufer(1:(N_box_Vx_i_top+N_box_Vx_i_top))  = i_vzdf_loc(N_box_Vx_i_low:N_box_Vx_i_top, loc)
!          ibufer2(1:(N_box_Vx_i_top+N_box_Vx_i_top)) = 0
!          CALL MPI_REDUCE(ibufer, ibufer2, (N_box_Vx_i_top+N_box_Vx_i_top), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
!          CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO

   ELSE                                                                                     ! server >>>

! receive data from clients
! N_in_loc
      ibufer  = 0
      ibufer2 = 0
      CALL MPI_REDUCE(ibufer2, ibufer, N_of_all_vdf_locs, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      N_in_loc(1:N_of_all_vdf_locs) = ibufer(1:N_of_all_vdf_locs)
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
! distributions . . .
      DO loc = 1, N_of_all_vdf_locs          ! cycle over locations
! i_vxdf_loc
         ibufer  = 0
         ibufer2 = 0
         CALL MPI_REDUCE(ibufer2, ibufer, (N_box_Vx_i_top+N_box_Vx_i_top), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         i_vxdf_loc(N_box_Vx_i_low:N_box_Vx_i_top, loc) = ibufer(1:(N_box_Vx_i_top+N_box_Vx_i_top))
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
! ! i_vydf_loc
!          ibufer  = 0
!          ibufer2 = 0
!          CALL MPI_REDUCE(ibufer2, ibufer, (N_box_Vx_i_top+N_box_Vx_i_top), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
!          i_vydf_loc(N_box_Vx_i_low:N_box_Vx_i_top, loc) = ibufer(1:(N_box_Vx_i_top+N_box_Vx_i_top))
!          CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
! ! i_vzdf_loc
!          ibufer  = 0
!          ibufer2 = 0
!          CALL MPI_REDUCE(ibufer2, ibufer, (N_box_Vx_i_top+N_box_Vx_i_top), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
!          i_vzdf_loc(N_box_Vx_i_low:N_box_Vx_i_top, loc) = ibufer(1:(N_box_Vx_i_top+N_box_Vx_i_top))
!          CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO

! IVDF(VX)
      ivxdf_filename = '_TTTT_ivxdf.dat'
      ivxdf_filename(2:5) = snapnumber_txt
      OPEN (99, FILE = ivxdf_filename)
      WRITE (99, '("# col   1 is Vx [in units of V_th_e = ",e12.5," m/s]")') V_Te_ms
      WRITE (99, '("#----------")')
      DO loc = 1, N_of_all_vdf_locs
         WRITE (99, '("# col ",i3," is the ion VDF(Vx) in location ",i2," with ",i8," macroparticles")') loc+1, loc, N_in_loc(loc)
      END DO
      WRITE (99, '("#----------")')
      DO indx_x = N_box_Vx_i_low, N_box_Vx_i_top
         temp_arr(1:N_of_all_vdf_locs) = DBLE(i_vxdf_loc(indx_x, 1:N_of_all_vdf_locs)) * DBLE(N_box_vel) * SQRT(Ms(2))
         WRITE (99, '(2x,f10.5,99(1x,e12.5))') ivx_mid_of_box(indx_x), temp_arr(1:N_of_all_vdf_locs)
      END DO
      CLOSE (99, STATUS = 'KEEP')
   END IF

   i_vxdf_loc = 0
   ! i_vydf_loc = 0
   ! i_vzdf_loc = 0

   IF (ALLOCATED(ibufer)) THEN
      DEALLOCATE(ibufer, STAT = DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : SNAP_LOCAL_ION_VDFS : Error in DEALLOCATE ibufer !!!")', Rank_of_process
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(ibufer2)) THEN
      DEALLOCATE(ibufer2, STAT = DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : SNAP_LOCAL_ION_VDFS : Error in DEALLOCATE ibufer2 !!!")', Rank_of_process
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

END SUBROUTINE SNAP_LOCAL_ION_VDFS

!------------------------------------------
! produces the phase plane X-VX for ions
SUBROUTINE SNAP_ION_PHASE_PLANE

   use mpi
   USE ParallelOperationValues
   USE Snapshots
   USE CurrentProblemValues

   IMPLICIT NONE

!  INCLUDE 'mpif.h'

   INTEGER k               ! particle index
   INTEGER n

   CHARACTER(13) ipp_filename

   REAL(8), ALLOCATABLE :: dbufer(:)
   INTEGER bufer_len, transm_len
   INTEGER ALLOC_ERR, DEALLOC_ERR, ierr
   INTEGER stattus(MPI_STATUS_SIZE)
!  LOGICAL flag

   IF (Flag_pp_creation.EQ.0) RETURN

   IF (Rank_of_process.GT.0) THEN                                                                             ! client >>>

      bufer_len = 4*((N_part(2) - 1) / (N_to_skip_ion + 1) + 10)     ! define the bufer length

      IF (.NOT.ALLOCATED(dbufer)) THEN
         ALLOCATE(dbufer(1:bufer_len), STAT=ALLOC_ERR)
         IF(ALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : SNAP_ION_PHASE_PLANE : Error in ALLOCATE dbufer !!!")', Rank_of_process
            PRINT '(2x,"The program will be terminated now :(")'
            STOP
         END IF
      END IF

      n = 1
      DO k = 1, N_part(2), N_to_skip_ion + 1                          ! fill the buffer with the particles data
         dbufer(n)   =  X_of_spec(2)%part(k)
         dbufer(n+1) = VX_of_spec(2)%part(k)
         dbufer(n+2) = VY_of_spec(2)%part(k)
         dbufer(n+3) = VZ_of_spec(2)%part(k)
         n = n + 4
      END DO
      transm_len = n - 1      ! now transm_len is the length of the transmission

      IF (transm_len.GT.bufer_len) THEN
         PRINT '(2x,"Process ",i3," : Error in SNAP_ION_PHASE_PLANES !!!")', Rank_of_process
         PRINT '(2x,"Index in the bufer ",i8," exceeds the threshold value ",i8)', transm_len, bufer_len
         PRINT '(2x,"The program will be terminated now :(")'
         STOP
      END IF

! transmit data
! bufer_len
      CALL MPI_SEND(transm_len, 1, MPI_INTEGER, 0, 21, MPI_COMM_WORLD, ierr)
! dbufer
      CALL MPI_SEND(dbufer, transm_len, MPI_DOUBLE_PRECISION, 0, 22, MPI_COMM_WORLD, ierr)

!     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

      IF (ALLOCATED(dbufer)) THEN
         DEALLOCATE(dbufer, STAT = DEALLOC_ERR)
         IF(DEALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : SNAP_ION_PHASE_PLANE : Error in DEALLOCATE dbufer !!!")', Rank_of_process
            PRINT *, 'The program will be terminated now :('
            STOP
         END IF
      END IF

   ELSE                                                                           ! server >>>>

      PRINT '(2x,"Process ",i3," : Creating ion phase plane...")', Rank_of_process

! Produce new filename
      ipp_filename =  '_TTTT_ipp.dat'
      ipp_filename(2:5) = snapnumber_txt
      OPEN (41, FILE = ipp_filename)

      DO k = 1, N_of_processes - 1

         CALL MPI_PROBE(k, 21, MPI_COMM_WORLD, stattus, ierr)
         CALL MPI_RECV(bufer_len, 1, MPI_INTEGER, k, 21, MPI_COMM_WORLD, stattus, ierr)

         IF (.NOT.ALLOCATED(dbufer)) THEN
            ALLOCATE(dbufer(1:bufer_len), STAT=ALLOC_ERR)
            IF(ALLOC_ERR.NE.0)THEN
               PRINT '(2x,"Process ",i3," : SNAP_ION_PHASE_PLANE : Error in ALLOCATE dbufer !!!")', Rank_of_process
               PRINT '(2x,"The program will be terminated now :(")'
               STOP
            END IF
         END IF

         dbufer = 0.0_8

         CALL MPI_PROBE(k, 22, MPI_COMM_WORLD, stattus, ierr)
         CALL MPI_RECV(dbufer, bufer_len, MPI_DOUBLE_PRECISION, k, 22, MPI_COMM_WORLD, stattus, ierr)

         DO n = 1, bufer_len - 3, 4
            WRITE (41, '(4(2x,e12.5))') dbufer(n), dbufer(n+1), dbufer(n+2), dbufer(n+3)
         END DO

         IF (ALLOCATED(dbufer)) THEN
            DEALLOCATE(dbufer, STAT = DEALLOC_ERR)
            IF(DEALLOC_ERR.NE.0)THEN
               PRINT '(2x,"Process ",i3," : SNAP_ION_PHASE_PLANE : Error in DEALLOCATE dbufer !!!")', Rank_of_process
               PRINT *, 'The program will be terminated now :('
               STOP
            END IF
         END IF

      END DO

      CLOSE (41, STATUS = 'KEEP')

!     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

   END IF

END SUBROUTINE SNAP_ION_PHASE_PLANE

!---------------------------------------------------------------------------------------------------
! adds electron, which hit the wall, to the corresponding 2-d distribution function of the left wall
SUBROUTINE ADD_PRIMARY_E_TO_LEFT_DF(Vx, Vy, Vz)

   USE Snapshots
!  USE CurrentProblemValues, ONLY : N_box_Vx_e, N_box_Vyz_e

   IMPLICIT NONE

   REAL(8) Vx       ! x-component of electron velocity  (dim-less)
   REAL(8) Vy       ! y-component of electron velocity  (dim-less)
   REAL(8) Vz       ! z-component of electron velocity  (dim-less)

   INTEGER index_norm   ! index of box for x-velocity
   INTEGER index_par    ! index of box for parallel velocity

! quit if the flag is not set for saving
   IF (.NOT.Accumulate_wall_df) RETURN

   index_norm = INT(ABS(Vx))
   index_par  = INT(SQRT(Vy*Vy+Vz*Vz))

   IF (index_norm.GT.N_box_Vx_e) RETURN
   IF (index_par.GT.N_box_Vyz_e) RETURN

   ep_2vdf_lw(index_norm, index_par) = ep_2vdf_lw(index_norm, index_par) + 1

END SUBROUTINE ADD_PRIMARY_E_TO_LEFT_DF

!----------------------------------------------------------------------------------------------------
! adds electron, which hit the wall, to the corresponding 2-d distribution function of the right wall
SUBROUTINE ADD_PRIMARY_E_TO_RIGHT_DF(Vx, Vy, Vz)

   USE Snapshots
!  USE CurrentProblemValues, ONLY : N_box_Vx_e, N_box_Vyz_e

   IMPLICIT NONE

   REAL(8) Vx       ! x-component of electron velocity  (dim-less)
   REAL(8) Vy       ! y-component of electron velocity  (dim-less)
   REAL(8) Vz       ! z-component of electron velocity  (dim-less)

   INTEGER index_norm   ! index of box for x-velocity
   INTEGER index_par    ! index of box for parallel velocity

! quit if the flag is not set for saving
   IF (.NOT.Accumulate_wall_df) RETURN

   index_norm = INT(ABS(Vx))
   index_par  = INT(SQRT(Vy*Vy+Vz*Vz))

   IF (index_norm.GT.N_box_Vx_e) RETURN
   IF (index_par.GT.N_box_Vyz_e) RETURN

   ep_2vdf_rw(index_norm, index_par) = ep_2vdf_rw(index_norm, index_par) + 1

END SUBROUTINE ADD_PRIMARY_E_TO_RIGHT_DF

!----------------------------------------------------------------------------------------------------------------
! adds electron, which was emitted from the wall, to the corresponding 2-d distribution function of the left wall
SUBROUTINE ADD_EMITTED_E_TO_LEFT_DF(Vx, Vy, Vz)

   USE Snapshots
!  USE CurrentProblemValues, ONLY : N_box_Vx_e, N_box_Vyz_e

   IMPLICIT NONE

   REAL(8) Vx       ! x-component of electron velocity  (dim-less)
   REAL(8) Vy       ! y-component of electron velocity  (dim-less)
   REAL(8) Vz       ! z-component of electron velocity  (dim-less)

   INTEGER index_norm   ! index of box for x-velocity
   INTEGER index_par    ! index of box for parallel velocity

! quit if the flag is not set for saving
   IF (.NOT.Accumulate_wall_df) RETURN

   index_norm = INT(ABS(Vx))
   index_par  = INT(SQRT(Vy*Vy+Vz*Vz))

   IF (index_norm.GT.N_box_Vx_e) RETURN
   IF (index_par.GT.N_box_Vyz_e) RETURN

   es_2vdf_lw(index_norm, index_par) = es_2vdf_lw(index_norm, index_par) + 1

END SUBROUTINE ADD_EMITTED_E_TO_LEFT_DF

!-----------------------------------------------------------------------------------------------------------------
! adds electron, which was emitted from the wall, to the corresponding 2-d distribution function of the right wall
SUBROUTINE ADD_EMITTED_E_TO_RIGHT_DF(Vx, Vy, Vz)

   USE Snapshots
!  USE CurrentProblemValues, ONLY : N_box_Vx_e, N_box_Vyz_e

   IMPLICIT NONE

   REAL(8) Vx       ! x-component of electron velocity  (dim-less)
   REAL(8) Vy       ! y-component of electron velocity  (dim-less)
   REAL(8) Vz       ! z-component of electron velocity  (dim-less)

   INTEGER index_norm   ! index of box for x-velocity
   INTEGER index_par    ! index of box for parallel velocity

! quit if the flag is not set for saving
   IF (.NOT.Accumulate_wall_df) RETURN

   index_norm = INT(ABS(Vx))
   index_par  = INT(SQRT(Vy*Vy+Vz*Vz))

   IF (index_norm.GT.N_box_Vx_e) RETURN
   IF (index_par.GT.N_box_Vyz_e) RETURN

   es_2vdf_rw(index_norm, index_par) = es_2vdf_rw(index_norm, index_par) + 1

END SUBROUTINE ADD_EMITTED_E_TO_RIGHT_DF

!-------------------------------------
!
SUBROUTINE CREATE_DF_ARRAYS

   USE CurrentProblemValues, ONLY : N_box_vel, V_Te_ms, Ms, N_spec
   USE Snapshots
   USE ParallelOperationValues

   IMPLICIT NONE
   INTEGER ALLOC_ERR, i
!  N_box_Vx_e_low  = -N_box_Vx_e
!  N_box_Vx_e_top  =  N_box_Vx_e + 1
!  N_box_Vyz_e_low = -N_box_Vyz_e
!  N_box_Vyz_e_top =  N_box_Vyz_e + 1
!  N_box_Vx_i_low  = -N_box_Vx_i
!  N_box_Vx_i_top  =  N_box_Vx_i + 1

! create and store the middles of the velocity boxes for the velocity distributions (in units of electron thermal velocity !!!)
   ! Note: we do not allocate velocity mid of box arrays based on distribution function
   !       flags in order to ensure they exist for reference distribution function creation.
   !       This could of course be changed in the future, if you want to make reference
   !       distribution functions flagged as well. Cheers.
   ALLOCATE(evx_mid_of_box(N_box_Vx_e_low:N_box_Vx_e_top), STAT=ALLOC_ERR)
   IF(ALLOC_ERR.NE.0)THEN
      PRINT *, 'Error in ALLOCATE evx_mid_of_box !!!'
      PRINT *, 'The program will be terminated now :('
      STOP
   END IF

   ALLOCATE(evyz_mid_of_box(N_box_Vyz_e_low:N_box_Vyz_e_top), STAT=ALLOC_ERR)
   IF(ALLOC_ERR.NE.0)THEN
      PRINT *, 'Error in ALLOCATE evyz_mid_of_box !!!'
      PRINT *, 'The program will be terminated now :('
      STOP
   END IF

   DO i = N_box_Vx_e_low, N_box_Vx_e_top
      evx_mid_of_box(i) = (DBLE(i) - 0.5_8) / DBLE(N_box_vel)
   END DO

   DO i = N_box_Vyz_e_low, N_box_Vyz_e_top
      evyz_mid_of_box(i) = (DBLE(i) - 0.5_8) / DBLE(N_box_vel)
   END DO

   ! create and store the middles of the velocity boxes for the velocity distributions (in units of electron thermal velocity !!!)
   if (flag_eedf) then
      ALLOCATE(eedf_mid_of_box_eV(1:N_E_bins), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in ALLOCATE eedf_mid_of_box_eV !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF

      DO i = 1, N_E_bins
         eedf_mid_of_box_eV(i) = (DBLE(i) - 0.5_8) * delta_Ee_eV
      END DO
   end if

!  IF (Rank_of_process.EQ.0) THEN
!! the two files below are necessary for processing of distributions, depending only on the absolute value of electron velocity
!     OPEN (99, FILE = 'dim_mid_evnorm.dat')
!     DO i = 1, N_box_Vx_e_top
!        WRITE (99, '(2x,i6,2x,f10.5,2x,e12.5)') i-1, evx_mid_of_box(i),  evx_mid_of_box(i) * V_Te_ms
!     END DO
!     CLOSE (99, STATUS = 'KEEP')
!
!     OPEN (99, FILE = 'dim_mid_evpar.dat')
!     DO i = 1, N_box_Vyz_e_top
!        WRITE (99, '(2x,i6,2x,f10.5,2x,e12.5)') i-1, evyz_mid_of_box(i),  evyz_mid_of_box(i) * V_Te_ms
!     END DO
!     CLOSE (99, STATUS = 'KEEP')
!  END IF

!----------
   ALLOCATE(ep_2vdf_lw(0:N_box_Vx_e,0:N_box_Vyz_e), STAT=ALLOC_ERR)
   IF(ALLOC_ERR.NE.0)THEN
      PRINT *, 'Error in ALLOCATE ep_2vdf_lw !!!'
      PRINT *, 'The program will be terminated now :('
      STOP
   END IF

   ALLOCATE(ep_2vdf_rw(0:N_box_Vx_e,0:N_box_Vyz_e), STAT=ALLOC_ERR)
   IF(ALLOC_ERR.NE.0)THEN
      PRINT *, 'Error in ALLOCATE ep_2vdf_rw !!!'
      PRINT *, 'The program will be terminated now :('
      STOP
   END IF

   ALLOCATE(es_2vdf_lw(0:N_box_Vx_e,0:N_box_Vyz_e), STAT=ALLOC_ERR)
   IF(ALLOC_ERR.NE.0)THEN
      PRINT *, 'Error in ALLOCATE es_2vdf_lw !!!'
      PRINT *, 'The program will be terminated now :('
      STOP
   END IF

   ALLOCATE(es_2vdf_rw(0:N_box_Vx_e,0:N_box_Vyz_e), STAT=ALLOC_ERR)
   IF(ALLOC_ERR.NE.0)THEN
      PRINT *, 'Error in ALLOCATE es_2vdf_rw !!!'
      PRINT *, 'The program will be terminated now :('
      STOP
   END IF

! clear the distribution function arrays
   ep_2vdf_lw = 0
   ep_2vdf_rw = 0
   es_2vdf_lw = 0
   es_2vdf_rw = 0

   IF (N_of_all_vdf_locs.EQ.0) RETURN    ! skip the rest if no local distributions are requested (N_of_breakpoints < 0 in the initialization file)

!  IF (Rank_of_process.EQ.0) THEN
!! the two files below are necessary only if the local distribution functions, using both positive and negative velocity, will be created
!     OPEN (99, FILE = 'dim_mid_evx.dat')
!     DO i = N_box_Vx_e_low, N_box_Vx_e_top
!        WRITE (99, '(2x,i6,2x,f10.5,2x,e12.5)') i, evx_mid_of_box(i),  evx_mid_of_box(i) * V_Te_ms
!     END DO
!     CLOSE (99, STATUS = 'KEEP')
!
!     OPEN (99, FILE = 'dim_mid_evyz.dat')
!     DO i = N_box_Vyz_e_low, N_box_Vyz_e_top
!        WRITE (99, '(2x,i6,2x,f10.5,2x,e12.5)') i, evyz_mid_of_box(i),  evyz_mid_of_box(i) * V_Te_ms
!     END DO
!     CLOSE (99, STATUS = 'KEEP')
!  END IF

! note, prevously there was only one 2d evdf array e_2vdf_loc(N_box_Vx_e_low:N_box_Vx_e_top,0:N_box_Vyz_e,1:N_of_all_vdf_locs)

!--------- common
!  ALLOCATE(e_2vxvydf_loc(N_box_Vx_e_low:N_box_Vx_e_top, N_box_Vyz_e_low:N_box_Vyz_e_top, 1:N_of_all_vdf_locs), STAT=ALLOC_ERR)
!  IF(ALLOC_ERR.NE.0)THEN
!     PRINT *, 'Error in ALLOCATE e_2vxvydf_loc !!!'
!     PRINT *, 'The program will be terminated now :('
!     STOP
!  END IF

!  ALLOCATE(e_2vxvzdf_loc(N_box_Vx_e_low:N_box_Vx_e_top, N_box_Vyz_e_low:N_box_Vyz_e_top, 1:N_of_all_vdf_locs), STAT=ALLOC_ERR)
!  IF(ALLOC_ERR.NE.0)THEN
!     PRINT *, 'Error in ALLOCATE e_2vxvzdf_loc !!!'
!     PRINT *, 'The program will be terminated now :('
!     STOP
!  END IF
   if (flag_evxdf) then
      ALLOCATE(e_vxdf_loc(N_box_Vx_e_low:N_box_Vx_e_top,1:N_of_all_vdf_locs), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in ALLOCATE e_vxdf_loc !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
      e_vxdf_loc = 0
   end if

   if (flag_evydf) then
      ALLOCATE(e_vydf_loc(N_box_Vyz_e_low:N_box_Vyz_e_top,1:N_of_all_vdf_locs), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in ALLOCATE e_vydf_loc !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
      e_vydf_loc = 0
   end if

   if (flag_evzdf) then
      ALLOCATE(e_vzdf_loc(N_box_Vyz_e_low:N_box_Vyz_e_top,1:N_of_all_vdf_locs), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in ALLOCATE e_vzdf_loc !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
      e_vzdf_loc = 0
   end if

   if (flag_eedf) then
      ALLOCATE(e_edf_loc(1:N_E_bins,1:N_of_all_edf_locs), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in ALLOCATE e_edf_loc !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
      e_edf_loc = 0
   end if

!--------- emitted from the left
!  ALLOCATE(ebl_2vxvydf_loc(N_box_Vx_e_low:N_box_Vx_e_top, N_box_Vyz_e_low:N_box_Vyz_e_top, 1:N_of_all_vdf_locs), STAT=ALLOC_ERR)
!  IF(ALLOC_ERR.NE.0)THEN
!     PRINT *, 'Error in ALLOCATE ebl_2vxvydf_loc !!!'
!     PRINT *, 'The program will be terminated now :('
!     STOP
!  END IF

!  ALLOCATE(ebl_2vxvzdf_loc(N_box_Vx_e_low:N_box_Vx_e_top, N_box_Vyz_e_low:N_box_Vyz_e_top, 1:N_of_all_vdf_locs), STAT=ALLOC_ERR)
!  IF(ALLOC_ERR.NE.0)THEN
!     PRINT *, 'Error in ALLOCATE ebl_2vxvzdf_loc !!!'
!     PRINT *, 'The program will be terminated now :('
!     STOP
!  END IF

   ALLOCATE(ebl_vxdf_loc(N_box_Vx_e_low:N_box_Vx_e_top,1:N_of_all_vdf_locs), STAT=ALLOC_ERR)
   IF(ALLOC_ERR.NE.0)THEN
      PRINT *, 'Error in ALLOCATE ebl_vxdf_loc !!!'
      PRINT *, 'The program will be terminated now :('
      STOP
   END IF

   ALLOCATE(ebl_vydf_loc(N_box_Vyz_e_low:N_box_Vyz_e_top,1:N_of_all_vdf_locs), STAT=ALLOC_ERR)
   IF(ALLOC_ERR.NE.0)THEN
      PRINT *, 'Error in ALLOCATE ebl_vydf_loc !!!'
      PRINT *, 'The program will be terminated now :('
      STOP
   END IF

   ALLOCATE(ebl_vzdf_loc(N_box_Vyz_e_low:N_box_Vyz_e_top,1:N_of_all_vdf_locs), STAT=ALLOC_ERR)
   IF(ALLOC_ERR.NE.0)THEN
      PRINT *, 'Error in ALLOCATE ebl_vzdf_loc !!!'
      PRINT *, 'The program will be terminated now :('
      STOP
   END IF

!--------- emitted from the right
!  ALLOCATE(ebr_2vxvydf_loc(N_box_Vx_e_low:N_box_Vx_e_top, N_box_Vyz_e_low:N_box_Vyz_e_top, 1:N_of_all_vdf_locs), STAT=ALLOC_ERR)
!  IF(ALLOC_ERR.NE.0)THEN
!     PRINT *, 'Error in ALLOCATE ebr_2vxvydf_loc !!!'
!     PRINT *, 'The program will be terminated now :('
!     STOP
!  END IF

!  ALLOCATE(ebr_2vxvzdf_loc(N_box_Vx_e_low:N_box_Vx_e_top, N_box_Vyz_e_low:N_box_Vyz_e_top, 1:N_of_all_vdf_locs), STAT=ALLOC_ERR)
!  IF(ALLOC_ERR.NE.0)THEN
!     PRINT *, 'Error in ALLOCATE ebr_2vxvzdf_loc !!!'
!     PRINT *, 'The program will be terminated now :('
!     STOP
!  END IF

   ALLOCATE(ebr_vxdf_loc(N_box_Vx_e_low:N_box_Vx_e_top,1:N_of_all_vdf_locs), STAT=ALLOC_ERR)
   IF(ALLOC_ERR.NE.0)THEN
      PRINT *, 'Error in ALLOCATE ebr_vxdf_loc !!!'
      PRINT *, 'The program will be terminated now :('
      STOP
   END IF

   ALLOCATE(ebr_vydf_loc(N_box_Vyz_e_low:N_box_Vyz_e_top,1:N_of_all_vdf_locs), STAT=ALLOC_ERR)
   IF(ALLOC_ERR.NE.0)THEN
      PRINT *, 'Error in ALLOCATE ebr_vydf_loc !!!'
      PRINT *, 'The program will be terminated now :('
      STOP
   END IF

   ALLOCATE(ebr_vzdf_loc(N_box_Vyz_e_low:N_box_Vyz_e_top,1:N_of_all_vdf_locs), STAT=ALLOC_ERR)
   IF(ALLOC_ERR.NE.0)THEN
      PRINT *, 'Error in ALLOCATE ebr_vzdf_loc !!!'
      PRINT *, 'The program will be terminated now :('
      STOP
   END IF

!  e_2vxvydf_loc = 0
!  e_2vxvzdf_loc = 0

!  ebl_2vxvydf_loc = 0
!  ebl_2vxvzdf_loc = 0
   ebl_vxdf_loc = 0
   ebl_vydf_loc = 0
   ebl_vzdf_loc = 0

!  ebr_2vxvydf_loc = 0
!  ebr_2vxvzdf_loc = 0
   ebr_vxdf_loc = 0
   ebr_vydf_loc = 0
   ebr_vzdf_loc = 0

!---------
   IF (N_spec.ge.2) THEN   ! if ions are accounted in simulation

      if (flag_ilwedf.or.flag_irwedf) then
         ALLOCATE(iedf_wall_mid_of_box_eV(1:N_E_bins), STAT=ALLOC_ERR)
         IF(ALLOC_ERR.NE.0)THEN
            PRINT *, 'Error in ALLOCATE iedf_wall_mid_of_box_eV !!!'
            PRINT *, 'The program will be terminated now :('
            STOP
         END IF
         DO i = 1, N_E_bins
            iedf_wall_mid_of_box_eV(i) = (DBLE(i) - 0.5_8) * delta_Ei_wall_eV
         END DO

         if (flag_irwedf) then
            ALLOCATE(irwedf(1 : N_E_bins), STAT=ALLOC_ERR)
            IF(ALLOC_ERR.NE.0)THEN
               PRINT *, 'Error in ALLOCATE irwedf !!!'
               PRINT *, 'The program will be terminated now :('
               STOP
            END IF
            irwedf = 0.0_8 !initialized here, cleared after each snapshot
         end if
   
         if (flag_ilwedf) then
            ALLOCATE(ilwedf(1 : N_E_bins), STAT=ALLOC_ERR)
            IF(ALLOC_ERR.NE.0)THEN
               PRINT *, 'Error in ALLOCATE ilwedf !!!'
               PRINT *, 'The program will be terminated now :('
               STOP
            END IF
            ilwedf = 0.0_8 !initialized here, cleared after each snapshot
         end if

      end if

      ALLOCATE(ivx_mid_of_box(N_box_Vx_i_low:N_box_Vx_i_top), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in ALLOCATE ivx_mid_of_box !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
      DO i = N_box_Vx_i_low, N_box_Vx_i_top
         ivx_mid_of_box(i) = (DBLE(i) - 0.5_8) / (DBLE(N_box_vel) * SQRT(Ms(2)))
      END DO

      if (flag_ivxdf) then
         ALLOCATE(i_vxdf_loc(N_box_Vx_i_low:N_box_Vx_i_top,1:N_of_all_vdf_locs), STAT=ALLOC_ERR)
         IF(ALLOC_ERR.NE.0)THEN
            PRINT *, 'Error in ALLOCATE i_vxdf_loc !!!'
            PRINT *, 'The program will be terminated now :('
            STOP
         END IF
         i_vxdf_loc = 0
      end if

      if (flag_iedf) then
         ALLOCATE(iedf_mid_of_box_eV(1:N_E_bins), STAT=ALLOC_ERR)
         IF(ALLOC_ERR.NE.0)THEN
            PRINT *, 'Error in ALLOCATE iedf_mid_of_box_eV !!!'
            PRINT *, 'The program will be terminated now :('
            STOP
         END IF
         DO i = 1, N_E_bins
            iedf_mid_of_box_eV(i) = (DBLE(i) - 0.5_8) * delta_Ei_eV
         END DO

         ALLOCATE(i_edf_loc(1:N_E_bins,1:N_of_all_edf_locs), STAT=ALLOC_ERR)
         IF(ALLOC_ERR.NE.0)THEN
            PRINT *, 'Error in ALLOCATE i_edf_loc !!!'
            PRINT *, 'The program will be terminated now :('
            STOP
         END IF
         i_edf_loc = 0
      end if

!     IF (Rank_of_process.EQ.0) THEN
!        OPEN (99, FILE = 'dim_mid_ivx.dat')
!        DO i = N_box_Vx_i_low, N_box_Vx_i_top
!           WRITE (99, '(2x,i6,2x,f10.5,2x,e12.5)') i, ivx_mid_of_box(i),  ivx_mid_of_box(i) * V_Te_ms
!        END DO
!        CLOSE (99, STATUS = 'KEEP')
!     END IF

      ! ALLOCATE(i_vydf_loc(N_box_Vx_i_low:N_box_Vx_i_top,1:N_of_all_vdf_locs), STAT=ALLOC_ERR)
      ! IF(ALLOC_ERR.NE.0)THEN
      !    PRINT *, 'Error in ALLOCATE i_vydf_loc !!!'
      !    PRINT *, 'The program will be terminated now :('
      !    STOP
      ! END IF
      ! i_vydf_loc = 0

      ! ALLOCATE(i_vzdf_loc(N_box_Vx_i_low:N_box_Vx_i_top,1:N_of_all_vdf_locs), STAT=ALLOC_ERR)
      ! IF(ALLOC_ERR.NE.0)THEN
      !    PRINT *, 'Error in ALLOCATE i_vzdf_loc !!!'
      !    PRINT *, 'The program will be terminated now :('
      !    STOP
      ! END IF
      ! i_vzdf_loc = 0

   END IF

END SUBROUTINE CREATE_DF_ARRAYS

!-------------------------------------
!
SUBROUTINE CREATE_STANDARD_DISTRIBUTIONS

   USE CurrentProblemValues, ONLY : N_box_vel, V_Te_ms, Ms, N_spec, T_e_eV, T_i_eV
   USE MCCollisions, ONLY : T_neutral_eV, Thresh_en_excit_eV, Thresh_en_ioniz_eV
   USE Snapshots

   IMPLICIT NONE

   INTEGER i       ! index of velocity box
   REAL(8) alfa    ! coefficient for species with different masses and temperatures
   REAL(8) alfa_i  ! is used for the electron component produced in the ionization
   REAL(8) df      ! distribution function value

! maxwellian distribution of electrons, T_e_eV
   alfa = 1.0_8
   OPEN (99, FILE = 'dim_std_maxwell_evdf.dat')
   DO i = N_box_Vx_e_low, N_box_Vx_e_top
      WRITE (99, '(2x,f9.4,2x,e12.5)') evx_mid_of_box(i), (1.0_8 / SQRT(3.1415926539_8)) * alfa * &
      & EXP(-(alfa * evx_mid_of_box(i))**2)
   END DO
   CLOSE (99, STATUS = 'KEEP')

   IF (N_spec.EQ.2) THEN

! maxwellian distribution of ions, T_i_eV
      alfa = SQRT((T_e_eV / T_i_eV) * Ms(2))
      OPEN (99, FILE = 'dim_std_maxwell_ivdf.dat')
      DO i = N_box_Vx_i_low, N_box_Vx_i_top
         WRITE (99, '(2x,f9.5,2x,e12.5)') ivx_mid_of_box(i), (1.0_8 / SQRT(3.1415926539_8)) * alfa * &
         & EXP(-(alfa * ivx_mid_of_box(i))**2)
      END DO
      CLOSE (99, STATUS = 'KEEP')

! maxwellian distribution of neutral gas, T_neutral_eV
      alfa = SQRT((T_e_eV / T_neutral_eV) * Ms(2))
      OPEN (99, FILE = 'dim_std_maxwell_nvdf.dat')
      DO i = N_box_Vx_i_low, N_box_Vx_i_top
         WRITE (99, '(2x,f9.5,2x,e12.5)') ivx_mid_of_box(i), (1.0_8 / SQRT(3.1415926539_8)) * alfa * &
         & EXP(-(alfa * ivx_mid_of_box(i))**2)
      END DO
      CLOSE (99, STATUS = 'KEEP')

   END IF

! 3d-isotropic distribution for monoenergetic electrons, T_e_eV
   OPEN (99, FILE = 'dim_std_iso3d_evdf.dat')
   DO i = N_box_Vx_e_low, N_box_Vx_e_top
      IF (ABS(evx_mid_of_box(i)).GT.1.0_8) THEN     ! 1.0_8 = VT(1)
         df = 0.0_8
      ELSE
         df = 0.5_8 / 1.0_8 ! 1.0_8 = VT(1)
      END IF
      WRITE (99, '(2x,f9.4,2x,e12.5)') evx_mid_of_box(i), df
   END DO
   CLOSE (99, STATUS = 'KEEP')

! 2d-isotropic distribution for monoenergetic electrons, T_e_eV
   OPEN (99, FILE = 'dim_std_iso2d_evdf.dat')
   DO i = N_box_Vx_e_low, N_box_Vx_e_top
      IF (ABS(evx_mid_of_box(i)).GE.1.0_8) THEN     ! 1.0_8 = VT(1)
         df = 0.0_8
      ELSE
         df = 1.0_8 / (3.1415926539_8 * SQRT(1.0_8 - evx_mid_of_box(i)**2))      !in denominator 1.0_8 = VT(1)**2
      END IF
      WRITE (99, '(2x,f9.4,2x,e12.5)') evx_mid_of_box(i), df
   END DO
   CLOSE (99, STATUS = 'KEEP')

! distribution for elastically scattered monoenergetic electron beam after scattering, over velocity parallel to the initial beam direction, energy T_e_eV
   OPEN (99, FILE = 'dim_std_scat_elast_evpardf.dat')
   alfa = T_e_eV
   DO i = N_box_Vyz_e_low, N_box_Vyz_e_top
      IF (ABS(evyz_mid_of_box(i)).GT.SQRT(alfa/T_e_eV)) THEN
         df = 0.0_8
      ELSE
         df = 1.0_8 / (LOG(1.0_8 + alfa) * (SQRT(alfa/T_e_eV) * (2.0_8 + alfa) / alfa - evyz_mid_of_box(i)))
      END IF
      WRITE (99, '(2x,f9.4,2x,e12.5)') evyz_mid_of_box(i), df
   END DO
   CLOSE (99, STATUS = 'KEEP')

! distribution for elastically scattered monoenergetic electron beam after scattering, over velocity normal to the initial beam direction, energy T_e_eV
   OPEN (99, FILE = 'dim_std_scat_elast_evnormdf.dat')
   alfa = T_e_eV
   DO i = N_box_Vyz_e_low, N_box_Vyz_e_top
      IF (ABS(evyz_mid_of_box(i)).GT.SQRT(alfa/T_e_eV)) THEN
         df = 0.0_8
      ELSE
         df = 1.0_8 / (LOG(1.0_8 + alfa) * SQRT(4.0_8 * (alfa/T_e_eV) * (1.0_8 + alfa) / alfa**2 + evyz_mid_of_box(i)**2))
      END IF
      WRITE (99, '(2x,f9.4,2x,e12.5)') evyz_mid_of_box(i), df
   END DO
   CLOSE (99, STATUS = 'KEEP')

! distribution for monoenergetic electron beam after scattering via the excitation, over velocity parallel to the initial beam direction,
! energy T_e_eV - Thresh_en_excit_eV
!  OPEN (99, FILE = 'dim_std_scat_excit_evpardf.dat')
!  alfa = T_e_eV - Thresh_en_excit_eV
!  DO i = N_box_Vyz_e_low, N_box_Vyz_e_top
!     IF (ABS(evyz_mid_of_box(i)).GT.SQRT(alfa/T_e_eV)) THEN
!        df = 0.0_8
!     ELSE
!        df = 1.0_8 / (LOG(1.0_8 + alfa) * (SQRT(alfa/T_e_eV) * (2.0_8 + alfa) / alfa - evyz_mid_of_box(i)))
!     END IF
!     WRITE (99, '(2x,f9.4,2x,e12.5)') evyz_mid_of_box(i), df
!  END DO
!  CLOSE (99, STATUS = 'KEEP')

! distribution for monoenergetic electron beam after scattering via the excitation, over velocity normal to the initial beam direction,
! energy T_e_eV - Thresh_en_excit_eV
!  OPEN (99, FILE = 'dim_std_scat_excit_evnormdf.dat')
!  alfa = T_e_eV - Thresh_en_excit_eV
!  DO i = N_box_Vyz_e_low, N_box_Vyz_e_top
!     IF (ABS(evyz_mid_of_box(i)).GT.SQRT(alfa/T_e_eV)) THEN
!        df = 0.0_8
!     ELSE
!        df = 1.0_8 / (LOG(1.0_8 + alfa) * SQRT(4.0_8 * (alfa/T_e_eV) * (1.0_8 + alfa) / alfa**2 + evyz_mid_of_box(i)**2))
!     END IF
!     WRITE (99, '(2x,f9.4,2x,e12.5)') evyz_mid_of_box(i), df
!  END DO
!  CLOSE (99, STATUS = 'KEEP')

! distribution for monoenergetic electron beam after scattering via the ionization, over velocity parallel to the initial beam direction,
! ASSUME that T_e_eV = 40eV, ionization yield is 100%,
! includes two components:  (i) scattered electrons with energy 5 eV and
!                          (ii) ejected electrons with energy T_e_eV - Thresh_en_excit_eV - 5 eV
   OPEN (99, FILE = 'dim_std_scat_ion40_evpardf.dat')
   alfa   = 5.0_8                                ! energy of incident electron after scattering
   alfa_i = 40.0_8 - Thresh_en_ioniz_eV - alfa   ! energy of produced electron
   DO i = N_box_Vyz_e_low, N_box_Vyz_e_top
! for scattered
      IF (ABS(evyz_mid_of_box(i)).GT.SQRT(alfa/40.0_8)) THEN
         df = 0.0_8
      ELSE
         df = 1.0_8 / (LOG(1.0_8 + alfa) * (SQRT(alfa/40.0_8) * (2.0_8 + alfa) / alfa - evyz_mid_of_box(i)))
      END IF
! for produced
      IF (ABS(evyz_mid_of_box(i)).GT.SQRT(alfa_i/40.0_8)) THEN
         df = df + 0.0_8
      ELSE
         df = df + 1.0_8 / (LOG(1.0_8 + alfa_i) * (SQRT(alfa_i/40.0_8) * (2.0_8 + alfa_i) / alfa_i - evyz_mid_of_box(i)))
      END IF
! conserve the normalization
      df = 0.5_8 * df
      WRITE (99, '(2x,f9.4,2x,e12.5)') evyz_mid_of_box(i), df
   END DO
   CLOSE (99, STATUS = 'KEEP')

! distribution for monoenergetic electron beam after scattering via the ionization, over velocity normal to the initial beam direction,
! ASSUME that T_e_eV = 40eV, ionization yield is 100%,
! includes two components:  (i) scattered electrons with energy 5 eV and
!                          (ii) ejected electrons with energy T_e_eV - Thresh_en_excit_eV - 5 eV
   OPEN (99, FILE = 'dim_std_scat_ion40_evnormdf.dat')
   alfa   = 5.0_8                                ! energy of incident electron after scattering
   alfa_i = 40.0_8 - Thresh_en_ioniz_eV - alfa   ! energy of produced electron
   DO i = N_box_Vyz_e_low, N_box_Vyz_e_top
! for scattered
      IF (ABS(evyz_mid_of_box(i)).GT.SQRT(alfa/40.0_8)) THEN
         df = 0.0_8
      ELSE
         df = 1.0_8 / (LOG(1.0_8 + alfa) * SQRT(4.0_8 * (alfa/T_e_eV) * (1.0_8 + alfa) / alfa**2 + evyz_mid_of_box(i)**2))
      END IF
! for produced
      IF (ABS(evyz_mid_of_box(i)).GT.SQRT(alfa_i/40.0_8)) THEN
         df = df + 0.0_8
      ELSE
         df = df + 1.0_8 / (LOG(1.0_8 + alfa_i) * SQRT(4.0_8 * (alfa_i/T_e_eV) * (1.0_8 + alfa_i) / alfa_i**2 + &
         & evyz_mid_of_box(i)**2))
      END IF
! conserve the normalization
      df = 0.5_8 * df
      WRITE (99, '(2x,f9.4,2x,e12.5)') evyz_mid_of_box(i), df
   END DO
   CLOSE (99, STATUS = 'KEEP')

END SUBROUTINE CREATE_STANDARD_DISTRIBUTIONS


!--------------------

SUBROUTINE FINISH_SNAPSHOTS

   USE Snapshots
   IMPLICIT NONE

   INTEGER DEALLOC_ERR

   IF (ALLOCATED(Tcntr_snapshot)) THEN
      DEALLOCATE(Tcntr_snapshot, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE Tcntr_snapshot !!!'
      END IF
   END IF

   IF (ALLOCATED(Vdf_location_bnd)) THEN
      DEALLOCATE(Vdf_location_bnd, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE Vdf_location_bnd !!!'
      END IF
   END IF

   IF (ALLOCATED(Edf_location_bnd)) THEN
      DEALLOCATE(Edf_location_bnd, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE Vdf_location_bnd !!!'
      END IF
   END IF
!-----------
   IF (ALLOCATED(evx_mid_of_box)) THEN
      DEALLOCATE(evx_mid_of_box, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE evx_mid_of_box !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(evyz_mid_of_box)) THEN
      DEALLOCATE(evyz_mid_of_box, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE evyz_mid_of_box !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(eedf_mid_of_box_eV)) THEN
      DEALLOCATE(eedf_mid_of_box_eV, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE eedf_mid_of_box_eV !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(iedf_mid_of_box_eV)) THEN
      DEALLOCATE(iedf_mid_of_box_eV, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE iedf_mid_of_box_eV !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(iedf_wall_mid_of_box_eV)) THEN
      DEALLOCATE(iedf_wall_mid_of_box_eV, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE iedf_wall_mid_of_box_eV !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF
!----------
   IF (ALLOCATED(ep_2vdf_lw)) THEN
      DEALLOCATE(ep_2vdf_lw, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE ep_2vdf_lw !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(ep_2vdf_rw)) THEN
      DEALLOCATE(ep_2vdf_rw, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE ep_2vdf_rw !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(es_2vdf_lw)) THEN
      DEALLOCATE(es_2vdf_lw, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE es_2vdf_lw !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(es_2vdf_rw)) THEN
      DEALLOCATE(es_2vdf_rw, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE es_2vdf_rw !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF
!!---------
!  IF (ALLOCATED(e_2vdf_loc)) THEN
!     DEALLOCATE(e_2vdf_loc, STAT=DEALLOC_ERR)
!     IF(DEALLOC_ERR.NE.0)THEN
!        PRINT *, 'Error in DEALLOCATE e_2vdf_loc !!!'
!        PRINT *, 'The program will be terminated now :('
!        STOP
!     END IF
!  END IF
!--------- common
   IF (ALLOCATED(e_vxdf_loc)) THEN
      DEALLOCATE(e_vxdf_loc, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE e_vxdf_loc !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(e_vydf_loc)) THEN
      DEALLOCATE(e_vydf_loc, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE e_vydf_loc !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(e_vzdf_loc)) THEN
      DEALLOCATE(e_vzdf_loc, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE e_vzdf_loc !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF
!--------- emitted from the left
   IF (ALLOCATED(ebl_vxdf_loc)) THEN
      DEALLOCATE(ebl_vxdf_loc, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE ebl_vxdf_loc !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(ebl_vydf_loc)) THEN
      DEALLOCATE(ebl_vydf_loc, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE ebl_vydf_loc !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(ebl_vzdf_loc)) THEN
      DEALLOCATE(ebl_vzdf_loc, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE ebl_vzdf_loc !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF
!--------- emitted from the right
   IF (ALLOCATED(ebr_vxdf_loc)) THEN
      DEALLOCATE(ebr_vxdf_loc, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE ebr_vxdf_loc !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(ebr_vydf_loc)) THEN
      DEALLOCATE(ebr_vydf_loc, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE ebr_vydf_loc !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(ebr_vzdf_loc)) THEN
      DEALLOCATE(ebr_vzdf_loc, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE ebr_vzdf_loc !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF
!-----------
   IF (ALLOCATED(irwedf)) THEN
      DEALLOCATE(irwedf, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE irwedf !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF
   IF (ALLOCATED(ilwedf)) THEN
      DEALLOCATE(ilwedf, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE ilwedf !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF


   IF (ALLOCATED(ivx_mid_of_box)) THEN
      DEALLOCATE(ivx_mid_of_box, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE ivx_mid_of_box !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(i_vxdf_loc)) THEN
      DEALLOCATE(i_vxdf_loc, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE i_vxdf_loc !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(e_edf_loc)) THEN
      DEALLOCATE(e_edf_loc, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE e_edf_loc !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(i_edf_loc)) THEN
      DEALLOCATE(i_edf_loc, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT *, 'Error in DEALLOCATE i_edf_loc !!!'
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

END SUBROUTINE FINISH_SNAPSHOTS

!------------------------------------------------------

SUBROUTINE SNAP_ION_EDF_RW_LW
!*** energy distribution of ions impacting the right wall 05/19/14
   use mpi
   USE ParallelOperationValues
   USE Snapshots
   USE CurrentProblemValues
   USE Diagnostics

   IMPLICIT NONE

!  INCLUDE 'mpif.h'

   INTEGER index_enr          ! index of velocity box

   CHARACTER(16) irwedf_filename, ilwedf_filename

   real(8), ALLOCATABLE :: rbufer(:)
   real(8), ALLOCATABLE :: rbufer2(:)
   INTEGER ALLOC_ERR, DEALLOC_ERR, ierr

   IF (.NOT.ALLOCATED(rbufer)) THEN
      ALLOCATE(rbufer(1:N_E_bins), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : SNAP_LOCAL_ION_VDFS : Error in ALLOCATE rbufer !!!")', Rank_of_process
         PRINT '(2x,"The program will be terminated now :(")'
         STOP
      END IF
   END IF

   IF (.NOT.ALLOCATED(rbufer2)) THEN
      ALLOCATE(rbufer2(1:N_E_bins), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : SNAP_LOCAL_ION_VDFS : Error in ALLOCATE rbufer2 !!!")', Rank_of_process
         PRINT '(2x,"The program will be terminated now :(")'
         STOP
      END IF
   END IF

   IF (Rank_of_process.GT.0) THEN  !transmit data to server
      if (flag_irwedf) then
         rbufer(1:N_E_bins)  = irwedf(1:N_E_bins)
         rbufer2(1:N_E_bins) = 0.0_8
         CALL MPI_REDUCE(rbufer, rbufer2, N_E_bins, MPI_double_precision, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      end if
      if (flag_ilwedf) then
         rbufer(1:N_E_bins)  = ilwedf(1:N_E_bins)
         rbufer2(1:N_E_bins) = 0.0_8
         CALL MPI_REDUCE(rbufer, rbufer2, N_E_bins, MPI_double_precision, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      end if
   ELSE  !receive data from clients
      if (flag_irwedf) then
         rbufer  = 0.0_8
         rbufer2 = 0.0_8
         CALL MPI_REDUCE(rbufer2, rbufer, N_E_bins, MPI_double_precision, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         irwedf(1:N_E_bins) = rbufer(1:N_E_bins)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      end if
      if (flag_ilwedf) then
         rbufer  = 0.0_8
         rbufer2 = 0.0_8
         CALL MPI_REDUCE(rbufer2, rbufer, N_E_bins, MPI_double_precision, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         ilwedf(1:N_E_bins) = rbufer(1:N_E_bins)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      end if
   END IF

   ! Print to the right wall df file for the current snapshot
   if (flag_irwedf) then
      irwedf_filename = '_TTTT_irwedf.dat'
      irwedf_filename(2:5) = snapnumber_txt
      OPEN (99, FILE = irwedf_filename)
      DO index_enr = 1, N_E_bins
         WRITE (99, '(2(1x,e12.5))') iedf_wall_mid_of_box_eV(index_enr), irwedf(index_enr)
      END DO
      CLOSE (99, STATUS = 'KEEP')
      irwedf = 0.0_8 !!!! cleared to accumulate again for the next snapshot, if any
   end if

   ! Print to the left wall df file for the current snapshot
   if (flag_ilwedf) then
      ilwedf_filename = '_TTTT_ilwedf.dat'
      ilwedf_filename(2:5) = snapnumber_txt
      OPEN (98, FILE = ilwedf_filename)
      DO index_enr = 1, N_E_bins
         WRITE (98, '(2(1x,e12.5))') iedf_wall_mid_of_box_eV(index_enr), ilwedf(index_enr)
      END DO
      CLOSE (98, STATUS = 'KEEP')
      ilwedf = 0.0_8 !!!! cleared to accumulate again for the next snapshot, if any
   end if

   IF (ALLOCATED(rbufer)) THEN
      DEALLOCATE(rbufer, STAT = DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : SNAP_ION_EDF_RW_LW: Error in DEALLOCATE rbufer !!!")', Rank_of_process
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(rbufer2)) THEN
      DEALLOCATE(rbufer2, STAT = DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : SNAP_ION_EDF_RW_LW: Error in DEALLOCATE rbufer2 !!!")', Rank_of_process
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

END SUBROUTINE SNAP_ION_EDF_RW_LW

!------------------------------------------------------
! produces the local ion energy distribution functions
SUBROUTINE SNAP_LOCAL_ION_EDFS

   use mpi
   USE ParallelOperationValues
   USE Snapshots
   USE CurrentProblemValues
   USE Diagnostics, ONLY : Factor_energy_eV

   IMPLICIT NONE

   !  INCLUDE 'mpif.h'

   INTEGER k               ! particle index
   INTEGER index_enr       ! index of energy box
   INTEGER loc             ! index of current location
   INTEGER cell            ! index of cell where the particle is

   REAL(8) v2              ! velocity magnitude squared (dim-less)
   REAL(8) energy_eV       ! particle energy in eV

   REAL(8) doublecell      ! particle position in double precision
   INTEGER endloc
   INTEGER N_all_temp

   CHARACTER(15) iedf_filename

   INTEGER N_in_loc(1:N_of_all_edf_locs)     ! number of particles inside the regions, where the distribution function is calculated
   ! these values are used for normalization
   ! As of now, we use the same locations for edf and vdf creation
   REAL(8) temp_arr(1:N_of_all_edf_locs)

   INTEGER, ALLOCATABLE :: ibufer(:)
   INTEGER, ALLOCATABLE :: ibufer2(:)
   INTEGER ALLOC_ERR, DEALLOC_ERR, ierr

   IF (N_of_all_edf_locs.LT.1) RETURN     ! quit if creation of local distributions was not requested

   IF (.NOT.ALLOCATED(ibufer)) THEN
      ALLOCATE(ibufer(1:N_E_bins), STAT=ALLOC_ERR)   ! Assume that 2 * N_E_bins is larger than the number of locations
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : SNAP_LOCAL_ION_EDFS : Error in ALLOCATE ibufer !!!")', Rank_of_process
         PRINT '(2x,"The program will be terminated now :(")'
         STOP
      END IF
   END IF

   IF (.NOT.ALLOCATED(ibufer2)) THEN
      ALLOCATE(ibufer2(1:N_E_bins), STAT=ALLOC_ERR)   ! Assume that 2 * N_E_bins is larger than the number of locations
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : SNAP_LOCAL_ION_EDFS : Error in ALLOCATE ibufer2 !!!")', Rank_of_process
         PRINT '(2x,"The program will be terminated now :(")'
         STOP
      END IF
   END IF

   IF (Rank_of_process.GT.0) THEN                         ! client >>>

      N_in_loc = 0

      ! calculate arrays of distribution functions
      DO k = 1, N_part(2)

         v2 = VX_of_spec(2)%part(k)**2 + VY_of_spec(2)%part(k)**2 + VZ_of_spec(2)%part(k)**2
         energy_eV = v2 * Ms(2) * Factor_energy_eV

         IF (energy_eV .le. Ei_max_eV) THEN
            index_enr = 1 + int(energy_eV/delta_Ei_eV)
            IF (index_enr .gt. N_E_bins) index_enr = N_E_bins
         END IF

         cell = INT(X_of_spec(2)%part(k))

         ! If cell greater than Edf_location_bnd(endloc), then the particle is outside the last location
         ! If this is the case, print a warning that includes the process number, and the particle number
         ! Then skip the particle and move on to the next one
         IF (cell.GE.Edf_location_bnd(N_of_all_edf_locs)) THEN
            PRINT '(2x,"Process ",i3," : SNAP_LOCAL_ION_EDFS : Warning: particle ",i8," is outside the last location !!!")', &
            & Rank_of_process, k
            PRINT '(2x,"The particle location is: cell = ",i8," !!!")', cell
            CYCLE
         END IF
         ! Note that this if statement above is not necessary if the particles are not allowed to leave the domain (which they are not)
         ! But, it seeemed like a compiler optimization error due to -O2 flags, so I added it just in case. Seems to work just fine now.
         !  - 9/25/23

         DO loc = 1, N_of_all_edf_locs
            IF (cell.LE.Edf_location_bnd(loc)) THEN

               i_edf_loc(index_enr, loc) = i_edf_loc(index_enr, loc) + 1
               N_in_loc(loc) = N_in_loc(loc) + 1

               EXIT
            END IF
         END DO

      END DO
      ! transmit data to the server
      ! N_in_loc
      ibufer(1:N_of_all_edf_locs)  = N_in_loc(1:N_of_all_edf_locs)
      ibufer2(1:N_of_all_edf_locs) = 0
      CALL MPI_REDUCE(ibufer, ibufer2, N_of_all_edf_locs, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

      ! distributions . . .
      DO loc = 1, N_of_all_edf_locs          ! cycle over locations
         ! i_edf_loc
         ibufer(1:N_E_bins)  = i_edf_loc(1:N_E_bins, loc)
         ibufer2(1:N_E_bins) = 0
         CALL MPI_REDUCE(ibufer, ibufer2, N_E_bins, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO

   ELSE                                                                                     ! server >>>

      ! receive data from clients
      ! N_in_loc
      ibufer  = 0
      ibufer2 = 0
      CALL MPI_REDUCE(ibufer2, ibufer, N_of_all_edf_locs, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      N_in_loc(1:N_of_all_edf_locs) = ibufer(1:N_of_all_edf_locs)
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      ! distributions . . .
      DO loc = 1, N_of_all_edf_locs          ! cycle over locations
         ! i_edf_loc
         ibufer  = 0
         ibufer2 = 0
         CALL MPI_REDUCE(ibufer2, ibufer, N_E_bins, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         i_edf_loc(1:N_E_bins, loc) = ibufer(1:N_E_bins)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO

      ! Write the EDFs
      iedf_filename = '_TTTT_iedf.dat'
      iedf_filename(2:5) = snapnumber_txt
      OPEN (99, FILE = iedf_filename)
      WRITE (99, '("# col   1 is the midpoint Energy [in units of eV]")')
      WRITE (99, '("#----------")')
      DO loc = 1, N_of_all_edf_locs
         WRITE (99, '("# col ",i3," is the ion EDF(Vx) in location ",i2," with ",i8," macroparticles")') loc+1, loc, N_in_loc(loc)
      END DO
      WRITE (99, '("#----------")')
      DO index_enr = 1, N_E_bins
         temp_arr(1:N_of_all_edf_locs) = DBLE(i_edf_loc(index_enr, 1:N_of_all_edf_locs))
         WRITE (99, '(1x,e13.6,99(1x,e12.5))') iedf_mid_of_box_eV(index_enr), temp_arr(1:N_of_all_edf_locs)
      END DO
      CLOSE (99, STATUS = 'KEEP')
   END IF

   i_edf_loc = 0 ! cleared to accumulate again for the next snapshot, if any

   IF (ALLOCATED(ibufer)) THEN
      DEALLOCATE(ibufer, STAT = DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : SNAP_LOCAL_ION_EDFS : Error in DEALLOCATE ibufer !!!")', Rank_of_process
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(ibufer2)) THEN
      DEALLOCATE(ibufer2, STAT = DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : SNAP_LOCAL_ION_EDFS : Error in DEALLOCATE ibufer2 !!!")', Rank_of_process
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

END SUBROUTINE SNAP_LOCAL_ION_EDFS

!------------------------------------------------------
! produces the local electron energy distribution functions
SUBROUTINE SNAP_LOCAL_ELECTRON_EDFS

   use mpi
   USE ParallelOperationValues
   USE Snapshots
   USE CurrentProblemValues
   USE Diagnostics, ONLY : Factor_energy_eV

   IMPLICIT NONE

   !  INCLUDE 'mpif.h'

   INTEGER k               ! particle index
   INTEGER index_enr       ! index of energy box
   INTEGER loc             ! index of current location
   INTEGER cell            ! index of cell where the particle is

   REAL(8) v2              ! velocity magnitude squared (dim-less)
   REAL(8) energy_eV       ! particle energy in eV

   CHARACTER(15) eedf_filename

   INTEGER N_in_loc(1:N_of_all_edf_locs)     ! number of particles inside the regions, where the distribution function is calculated
   ! these values are used for normalization
   ! As of now, we use the same locations for edf and vdf creation
   REAL(8) temp_arr(1:N_of_all_edf_locs)

   INTEGER, ALLOCATABLE :: ibufer(:)
   INTEGER, ALLOCATABLE :: ibufer2(:)
   INTEGER ALLOC_ERR, DEALLOC_ERR, ierr

   IF (N_of_all_edf_locs.LT.1) RETURN     ! quit if creation of local distributions was not requested

   IF (.NOT.ALLOCATED(ibufer)) THEN
      ALLOCATE(ibufer(1:N_E_bins), STAT=ALLOC_ERR)   ! Assume that 2 * N_E_bins is larger than the number of locations
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : SNAP_LOCAL_ELECTRON_EDFS : Error in ALLOCATE ibufer !!!")', Rank_of_process
         PRINT '(2x,"The program will be terminated now :(")'
         STOP
      END IF
   END IF

   IF (.NOT.ALLOCATED(ibufer2)) THEN
      ALLOCATE(ibufer2(1:N_E_bins), STAT=ALLOC_ERR)   ! Assume that 2 * N_E_bins is larger than the number of locations
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : SNAP_LOCAL_ELECTRON_EDFS : Error in ALLOCATE ibufer2 !!!")', Rank_of_process
         PRINT '(2x,"The program will be terminated now :(")'
         STOP
      END IF
   END IF

   IF (Rank_of_process.GT.0) THEN                         ! client >>>

      N_in_loc = 0

      ! calculate arrays of distribution functions
      DO k = 1, N_part(1)

         v2 = VX_of_spec(1)%part(k)**2 + VY_of_spec(1)%part(k)**2 + VZ_of_spec(1)%part(k)**2
         energy_eV = v2 * Ms(1) * Factor_energy_eV

         IF (energy_eV .le. Ee_max_eV) THEN
            index_enr = 1 + int(energy_eV/delta_Ee_eV)
            IF (index_enr .gt. N_E_bins) index_enr = N_E_bins
         END IF

         cell = INT(X_of_spec(1)%part(k))
         DO loc = 1, N_of_all_edf_locs
            IF (cell.LE.Edf_location_bnd(loc)) THEN

               e_edf_loc(index_enr, loc) = e_edf_loc(index_enr, loc) + 1
               N_in_loc(loc) = N_in_loc(loc) + 1

               EXIT
            END IF
         END DO

      END DO
      ! transmit data to the server
      ! N_in_loc
      ibufer(1:N_of_all_edf_locs)  = N_in_loc(1:N_of_all_edf_locs)
      ibufer2(1:N_of_all_edf_locs) = 0
      CALL MPI_REDUCE(ibufer, ibufer2, N_of_all_edf_locs, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

      ! distributions . . .
      DO loc = 1, N_of_all_edf_locs          ! cycle over locations
         ! e_edf_loc
         ibufer(1:N_E_bins)  = e_edf_loc(1:N_E_bins, loc)
         ibufer2(1:N_E_bins) = 0
         CALL MPI_REDUCE(ibufer, ibufer2, N_E_bins, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO

   ELSE                                                                                     ! server >>>

      ! receive data from clients
      ! N_in_loc
      ibufer  = 0
      ibufer2 = 0
      CALL MPI_REDUCE(ibufer2, ibufer, N_of_all_edf_locs, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      N_in_loc(1:N_of_all_edf_locs) = ibufer(1:N_of_all_edf_locs)
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      ! distributions . . .
      DO loc = 1, N_of_all_edf_locs          ! cycle over locations
         ! e_edf_loc
         ibufer  = 0
         ibufer2 = 0
         CALL MPI_REDUCE(ibufer2, ibufer, N_E_bins, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         e_edf_loc(1:N_E_bins, loc) = ibufer(1:N_E_bins)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      END DO

      ! Write the EDFs
      eedf_filename = '_TTTT_eedf.dat'
      eedf_filename(2:5) = snapnumber_txt
      OPEN (99, FILE = eedf_filename)
      WRITE (99, '("# col   1 is the midpoint Energy [in units of eV]")')
      WRITE (99, '("#----------")')
      DO loc = 1, N_of_all_edf_locs
         WRITE (99, '("# col ",i3," is the electron EDF(Vx) in location ",i2," with ",i8," macroparticles")') loc+1, loc, N_in_loc(loc)
      END DO
      WRITE (99, '("#----------")')
      DO index_enr = 1, N_E_bins
         temp_arr(1:N_of_all_edf_locs) = DBLE(e_edf_loc(index_enr, 1:N_of_all_edf_locs))
         WRITE (99, '(1x,e13.6,99(1x,e12.5))') eedf_mid_of_box_eV(index_enr), temp_arr(1:N_of_all_edf_locs)
      END DO
      CLOSE (99, STATUS = 'KEEP')
   END IF

   e_edf_loc = 0 ! cleared to accumulate again for the next snapshot, if any

   IF (ALLOCATED(ibufer)) THEN
      DEALLOCATE(ibufer, STAT = DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : SNAP_LOCAL_ELECTRON_EDFS : Error in DEALLOCATE ibufer !!!")', Rank_of_process
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(ibufer2)) THEN
      DEALLOCATE(ibufer2, STAT = DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : SNAP_LOCAL_ELECTRON_EDFS : Error in DEALLOCATE ibufer2 !!!")', Rank_of_process
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

END SUBROUTINE SNAP_LOCAL_ELECTRON_EDFS
