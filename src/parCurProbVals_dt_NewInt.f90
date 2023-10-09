!===================================================================================================
!
SUBROUTINE INITIATE_PARAMETERS

  use mpi
  USE ParallelOperationValues
  USE CurrentProblemValues
  USE Diagnostics
  USE Snapshots
  USE mt19937

  IMPLICIT NONE

!  INCLUDE 'mpif.h'

  LOGICAL exists
  INTEGER i, s
  INTEGER itmp

  INTEGER flag_rf                         ! flag for if rf is on or off at the left wall
  REAL(8) L_plasma_deb                    ! plasma length in debye lengths                  ! ####### ######### ######## was INTEGER ######## ###### ######

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  REAL RAN, R

  CHARACTER (77) buf
  CHARACTER(3)   procnumber_txt
  INTEGER proc_blnks

  INTEGER dummy_i
  REAL(8) dummy_r
  REAL(8) Energy_ebeam_eV

  REAL(8) W_cycl_min_s1, W_cycl_max_s1, W_cycl_temp_s1
  INTEGER i_Wc_min, i_Wc_max

! functions
  REAL(8) Bx_gauss
  REAL(8) By_gauss


  INQUIRE (FILE = 'ssc_initial.dat', EXIST = exists)
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  OPEN (9, FILE = 'ssc_initial.dat')

  IF(exists) THEN

     IF (Rank_of_process.EQ.0) THEN
        PRINT '(2x,"Process ",i3," : ssc_initial.dat is found. Reading the data file...")', Rank_of_process
     END IF

     READ (9, '(A77)') buf !============================ SYSTEM CONFIGURATION ===========================")')
     READ (9, '(A77)') buf !------------- BC: Left wall potential is ------------------------------------")')
     READ (9, '(A77)') buf !--            0 = fixed, or connected to RF if the next option is on       --")')
     READ (9, '(A77)') buf !--            1 = floating                                                 --")')
     READ (9, '(A77)') buf !--            2 = connected to external circuit                            --")')
     READ (9, '(A77)') buf !--            3 = Electrode current is specified                           --")')
     READ (9, '(A77)') buf !-------d----- 4 = RF potential applied directly or with dielectric layer ----")')
     READ (9, '(7x,i1)') BC_flag
     READ (9, '(A77)') buf !-------d----- Left wall is RF biased (1/0=Yes/No) ---------------------------")')
     READ (9, '(7x,i1)') flag_rf
     READ (9, '(A77)') buf !--dddddd.ddd- Resistance of the resistor in the external circuit (ohm) ------")')
     READ (9, '(2x,f10.3)') R_ext_ohm
     READ (9, '(A77)') buf !--dddddd.ddd- Electrode area (cm^2) -----------------------------------------")')
     READ (9, '(2x,f10.3)') S_electrode_cm2
     READ (9, '(A77)') buf !--dddddd.ddd- External z-electric field parallel to walls (V/m) -------------")')
     READ (9, '(2x,f10.3)') E_z_ext_Vm
     READ (9, '(A77)') buf !--dddddd.ddd- Left wall potential/battery voltage (V) -----------------------")')
     READ (9, '(2x,f10.3)') U_ext_V
     READ (9, '(A77)') buf !--#d.dddE#dd- Left wall RF frequency (Hz) -----------------------------------")')
     READ (9, '(2x,e10.3)') f_rf_Hz
     READ (9, '(A77)') buf !--dddddd.ddd- Left wall RF amplitude (V) ------------------------------------")')
     READ (9, '(2x,f10.3)') U_rf_V
     READ (9, '(A77)') buf !--#d.dddE#dd- Left wall RF start time (s) -----------------------------------")')
     READ (9, '(2x,e10.3)') t_start_s
     READ (9, '(A77)') buf !--dddddd.ddd- Plasma layer width (m, beam wavelengths if negative) ----------")')
     READ (9, '(2x,f11.4)') L_plasma_m
     READ (9, '(A77)') buf !--dddddd.ddd- Ion mass (a.m.u.) ---------------------------------------------")')
     READ (9, '(2x,f10.3)') M_i_amu
     READ (9, '(A77)') buf !--dddddd.ddd- Ion temperature (eV) ------------------------------------------")')
     READ (9, '(2x,f10.3)') T_i_eV
     READ (9, '(A77)') buf !======================= PARAMETERS DEFINING SCALE VALUES ====================")')
     READ (9, '(A77)') buf !--#d.dddE#dd- Plasma density (m^-3) -----------------------------------------")')
     READ (9, '(2x,e10.3)') N_plasma_m3
     READ (9, '(A77)') buf !--dddddd.ddd- Electron temperature (eV) -------------------------------------")')
     READ (9, '(2x,f10.3)') T_e_eV
     READ (9, '(A77)') buf !--dddddd----- Number of macroparticles per cell -----------------------------")')
     READ (9, '(2x,i6)') N_of_particles_cell
     READ (9, '(A77)') buf !--dddddd----- Number of cells per Debye length, micron_flag -----------------")')
     READ (9, '(A77)') buf !--            If micron_flag = 0:   delta_x found from scaling parameters ---")')
     READ (9, '(A77)') buf !--            If micron_flag > 0:   delta_x [mkm] = N_cells * 10^(n-1) ------")')
     READ (9, '(A77)') buf !--dddddd-#d-- If micron_flag < 0:   delta_x [mkm] = N_cells * 10^(n) --------")')
     READ (9, '(2x,i6,1x,i2)') N_of_cells_debye, micron_flag
     READ (9, '(A77)') buf !------------- Maximal expected velocity (in V_therm_e), picosec_flag --------")')
     READ (9, '(A77)') buf !--            If picosec_flag = 0:   delta_t found from scaling parameters --")')
     READ (9, '(A77)') buf !--            If picosec_flag > 0:   delta_t [ps] = Max_vel * 10^(n-1) ------")')
     READ (9, '(A77)') buf !--dddddd-#d-- If picosec_flag < 0:   delta_t [ps] = Max_vel * 10^(n) --------")')
     READ (9, '(2x,i6,1x,i2)') N_max_vel, picosec_flag
     READ (9, '(A77)') buf !--dddddd----- Number of velocity boxes per unit of V_therm ------------------")')
     READ (9, '(2x,i6)') N_box_vel
     READ (9, '(A77)') buf !============================ SIMULATION CONTROL =============================")')
     READ (9, '(A77)') buf !-------d----- Initial particle density profile (0/1=Uniform/Parabolic) ------")')
     READ (9, '(7x,i1)') Density_flag
     READ (9, '(A77)') buf !--dddddd.ddd- Duration of simulation (ns) -----------------------------------")')
     READ (9, '(2x,f10.3)') T_sim_ns   
     READ (9, '(A77)') buf !--dddddd----- Step for saving checkpoints (timesteps, skip if <=0) ----------")')
     READ (9, '(2x,i6)') SaveCheck_step
     READ (9, '(A77)') buf !--dddddd----- Seed for random numbers generator -----------------------------")')
     READ (9, '(2x,i6)') I_random_seed
     READ (9, '(A77)') buf !------------- Method of initialization --------------------------------------")')
     READ (9, '(A77)') buf !--            0 = start a new run, ordinary                                --")')
     READ (9, '(A77)') buf !--            1 = continue the old run, start at the last checkpoint       --")')
     READ (9, '(A77)') buf !-------d----- 2 = start a new run, take particles from the checkpoint -------")')
     READ (9, '(7x,i1)') Restore_from_checkpoint
           
  ELSE                                   ! if file does not exist

!======================= System configuration =====================
     BC_flag         = 0
     flag_rf         = 1
     R_ext_ohm       = 0.0_8
     S_electrode_cm2 = 400.0_8
     E_z_ext_Vm      = 10000.0_8
     U_ext_V         = 0.0_8      ! 
     f_rf_Hz         = 1.356d7
     U_rf_V          = 250.0_8
     t_start_s       = 0.0_8
     L_plasma_m      = 0.02_8
     M_i_amu         = 131.0_8    ! Xenon
     T_i_eV          = 1.0_8
!==================== Parameters defining scale values =============
     N_plasma_m3          = 1.0d17
     T_e_eV               = 40.0_8
     N_of_particles_cell  = 500
     N_of_cells_debye     = 16
     micron_flag          = 0
     N_max_vel            = 4
     picosec_flag         = 0
     N_box_vel            = 20
!======================= Simulation control ========================
     Density_flag            = 0
     T_sim_ns                = 1.0_8
     SaveCheck_step          = 200000
     I_random_seed           = 1234
     Restore_from_checkpoint = 0

     IF (Rank_of_process.EQ.0) THEN

        PRINT '(/2x,"Process ",i3," : File with the name ssc_initial.dat not found. Use the default settings ...")', Rank_of_process

        PRINT '(2x,"Process ",i3," : Create ssc_initial.dat file . . .")', Rank_of_process

        WRITE (9, '("============================ SYSTEM CONFIGURATION ===========================")')
        WRITE (9, '("------------- BC: Left wall potential is ------------------------------------")')
        WRITE (9, '("--            0 = fixed, or connected to RF if the next option is on       --")')
        WRITE (9, '("--            1 = floating                                                 --")')
        WRITE (9, '("--            2 = connected to external circuit                            --")')
        WRITE (9, '("--            3 = Electrode current is specified                           --")')
        WRITE (9, '("-------d----- 4 = RF potential applied directly or with dielectric layer ----")')
        WRITE (9, '(7x,i1)') BC_flag
        WRITE (9, '("-------d----- Left wall is RF biased (1/0=Yes/No) ---------------------------")')
        WRITE (9, '(7x,i1)') flag_rf
        WRITE (9, '("--dddddd.ddd- Resistor in the external circuit (ohm) ------------------------")')
        WRITE (9, '(2x,f10.3)') R_ext_ohm
        WRITE (9, '("--dddddd.ddd- Electrode area (cm^2) -----------------------------------------")')
        WRITE (9, '(2x,f10.3)') S_electrode_cm2
        WRITE (9, '("--dddddd.ddd- External z-electric field parallel to walls (V/m) -------------")')
        WRITE (9, '(2x,f10.3)') E_z_ext_Vm
        WRITE (9, '("--dddddd.ddd- Left wall potential/battery voltage (V) -----------------------")')
        WRITE (9, '(2x,f10.3)') U_ext_V
        WRITE (9, '("--#d.dddE#dd- Left wall RF frequency (Hz) -----------------------------------")')
        WRITE (9, '(2x,e10.3)') f_rf_Hz
        WRITE (9, '("--dddddd.ddd- Left wall RF amplitude (V) ------------------------------------")')
        WRITE (9, '(2x,f10.3)') U_rf_V
        WRITE (9, '("--#d.dddE#dd- Left wall RF start time (s) -----------------------------------")')
        WRITE (9, '(2x,e10.3)') t_start_s
        WRITE (9, '("--dddddd.ddd- Plasma layer width (m, beam wavelengths if negative) ----------")')
        WRITE (9, '(2x,f11.4)') L_plasma_m
        WRITE (9, '("--dddddd.ddd- Ion mass (a.m.u., m_e if negative) ----------------------------")')
        WRITE (9, '(2x,f10.3)') M_i_amu
        WRITE (9, '("--dddddd.ddd- Ion temperature (eV) ------------------------------------------")')
        WRITE (9, '(2x,f10.3)') T_i_eV
        WRITE (9, '("======================= PARAMETERS DEFINING SCALE VALUES ====================")')
        WRITE (9, '("--#d.dddE#dd- Plasma density (m^-3) -----------------------------------------")')
        WRITE (9, '(2x,e10.3)') N_plasma_m3
        WRITE (9, '("--dddddd.ddd- Electron temperature (eV) -------------------------------------")')
        WRITE (9, '(2x,f10.3)') T_e_eV
        WRITE (9, '("--dddddd----- Number of macroparticles per cell -----------------------------")')
        WRITE (9, '(2x,i6)') N_of_particles_cell
        WRITE (9, '("--dddddd----- Number of cells per Debye length, micron_flag -----------------")')
        WRITE (9, '("--            If micron_flag = 0:   delta_x found from scaling parameters ---")')
        WRITE (9, '("--            If micron_flag > 0:   delta_x [mkm] = N_cells * 10^(n-1) ------")')
        WRITE (9, '("--dddddd-#d-- If micron_flag < 0:   delta_x [mkm] = N_cells * 10^(n) --------")')
        WRITE (9, '(2x,i6,1x,i2)') N_of_cells_debye, micron_flag
        WRITE (9, '("------------- Maximal expected velocity (in V_therm_e), picosec_flag --------")')
        WRITE (9, '("--            If picosec_flag = 0:   delta_t found from scaling parameters --")')
        WRITE (9, '("--            If picosec_flag > 0:   delta_t [ps] = Max_vel * 10^(n-1) ------")')
        WRITE (9, '("--dddddd-#d-- If picosec_flag < 0:   delta_t [ps] = Max_vel * 10^(n) --------")')
        WRITE (9, '(2x,i6,1x,i2)') N_max_vel, picosec_flag
        WRITE (9, '("--dddddd----- Number of cells per Debye length ------------------------------")')
        WRITE (9, '(2x,i6)') N_of_cells_debye
        WRITE (9, '("----dddd----- Maximal expected velocity (in V_therm_e) ----------------------")')
        WRITE (9, '(4x,i4)') N_max_vel
        WRITE (9, '("--dddddd----- Number of velocity boxes per unit of V_therm ------------------")')
        WRITE (9, '(2x,i6)') N_box_vel
        WRITE (9, '("============================ SIMULATION CONTROL =============================")')
        WRITE (9, '("-------d----- Initial particle density profile (0/1=Uniform/Parabolic) ------")')
        WRITE (9, '(7x,i1)') Density_flag
        WRITE (9, '("--dddddd.ddd- Duration of simulation (ns) -----------------------------------")')
        WRITE (9, '(2x,f10.3)') T_sim_ns   
        WRITE (9, '("--dddddd----- Step for saving checkpoints (timesteps, skip if <=0) ----------")')
        WRITE (9, '(2x,i6)') SaveCheck_step
        WRITE (9, '("--dddddd----- Seed for random numbers generator -----------------------------")')
        WRITE (9, '(2x,i6)') I_random_seed
        WRITE (9, '("------------- Method of initialization --------------------------------------")')
        WRITE (9, '("--            0 = start a new run, ordinary                                --")')
        WRITE (9, '("--            1 = continue the old run, start at the last checkpoint       --")')
        WRITE (9, '("-------d----- 2 = start a new run, take particles from the checkpoint -------")')
        WRITE (9, '(7x,i1)') Restore_from_checkpoint

     END IF

  END IF

  CLOSE (9, STATUS = 'KEEP')

! fool-proofing
  IF ((BC_flag.EQ.2).AND.(R_ext_ohm.LE.0.0_8)) THEN
     IF (Rank_of_process.EQ.0) THEN
        PRINT '("INITIATE_PARAMETERS: ERROR: external contour boundary condition requires positive resistance")'
        PRINT '("the present value of R_ext_ohm is ",e12.5," ohm")', R_ext_ohm
        PRINT '("terminating the program")'
     END IF
     STOP
  END IF

  IF ((BC_flag.EQ.2).AND.(S_electrode_cm2.LE.0.0_8)) THEN
     IF (Rank_of_process.EQ.0) THEN
        PRINT '("INITIATE_PARAMETERS: ERROR: external contour boundary condition requires positive electrode area")'
        PRINT '("the present value of S_electrode_cm2 is ",e12.5," cm^2")', R_ext_ohm
        PRINT '("terminating the program")'
     END IF
     STOP
  END IF

  IF (flag_rf.EQ.0) THEN
     rf_on = .FALSE.
     ! Check if the BC condition 4 is used. If yes, then terminate the program
     ! NOTE: BC 4 with no rf may not break the code, but it has not been verified and was not intended to be used
     IF (BC_flag.EQ.4) THEN
       IF (Rank_of_process.EQ.0) THEN
         PRINT '("INITIATE_PARAMETERS: ERROR: BC condition 4 requres left wall rf to be on")'
         PRINT '("terminating the program")'
       END IF
       STOP
     END IF
  ELSE
     rf_on = .TRUE.
     IF (f_rf_Hz.LT.0.0_8) THEN
        IF (Rank_of_process.EQ.0) THEN
           PRINT '("INITIATE_PARAMETERS: ERROR: rf battery on the left wall requires positive rf frequency")'
           PRINT '("the present value of f_rf_Hz is ",e12.5," Hz")', f_rf_Hz
           PRINT '("terminating the program")'
        END IF
     STOP
     END IF
  END IF


  
  INQUIRE (FILE = 'ssc_anisotropy.dat', EXIST = exists)
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  OPEN (9, FILE = 'ssc_anisotropy.dat')

  IF (exists) THEN

     IF (Rank_of_process.EQ.0) THEN
        PRINT '(2x,"Process ",i3," : ssc_anisotropy.dat is found. Reading the data file...")', Rank_of_process
     END IF

     READ (9, '(A77)') buf !============== Parameters of anisotropic electron distribution ==============")')
     READ (9, '(A77)') buf !--dddddd.ddd- X-temperature (eV) --------------------------------------------")')
     READ (9, '(2x,f10.3)') Tx_e_eV
     READ (9, '(A77)') buf !--dddddd.ddd- Z-temperature (eV) --------------------------------------------")')
     READ (9, '(2x,f10.3)') Tz_e_eV
     READ (9, '(A77)') buf !------dd----- Maximal velocity for initial distribution (in V_therm_e) ------")')
     READ (9, '(6x,i2)') N_max_vel_distrib
     READ (9, '(A77)') buf !--#d.dddE#dd- Plasma density (m^-3), does not affect scaling ----------------")')
     READ (9, '(2x,e10.3)') N_distrib_m3
     READ (9, '(A77)') buf !=============== To calculate the drift velocity one must set: ===============")')
     READ (9, '(A77)') buf !--dddddd.ddd- E_z (V/m), not to be used in the equation of motion -----------")')
     READ (9, '(2x,f10.3)') Ez_drift_Vm
     READ (9, '(A77)') buf !--dddddd.ddd- B_x (Gauss), not to be used in the equation of motion ---------")')
     READ (9, '(2x,f10.3)') Bx_drift_Gs
           
  ELSE                                   ! if file does not exist

     Tx_e_eV           = T_e_eV 
     Tz_e_eV           = T_e_eV
     N_max_vel_distrib = 3
     N_distrib_m3      = N_plasma_m3
     Ez_drift_Vm       = 0.0_8
     Bx_drift_Gs       = 0.0_8

     IF (Rank_of_process.EQ.0) THEN

        PRINT '(/2x,"Process ",i3," : File with the name ssc_anisotropy.dat not found. Use the default settings ...")', Rank_of_process

        PRINT '(2x,"Process ",i3," : Create ssc_anisotropy.dat file . . .")', Rank_of_process

        WRITE (9, '("============== Parameters of anisotropic electron distribution ==============")')
        WRITE (9, '("--dddddd.ddd- X-temperature (eV) --------------------------------------------")')
        WRITE (9, '(2x,f10.3)') Tx_e_eV
        WRITE (9, '("--dddddd.ddd- Z-temperature (eV) --------------------------------------------")')
        WRITE (9, '(2x,f10.3)') Tz_e_eV
        WRITE (9, '("------dd----- Maximal velocity for initial distribution (in V_therm_e) ------")')
        WRITE (9, '(6x,i2)') N_max_vel_distrib
        WRITE (9, '("--#d.dddE#dd- Plasma density (m^-3), does not affect scaling ----------------")')
        WRITE (9, '(2x,e10.3)') N_distrib_m3
        WRITE (9, '("=============== To calculate the drift velocity one must set: ===============")')
        WRITE (9, '("--dddddd.ddd- E_z (V/m), not to be used in the equation of motion -----------")')
        WRITE (9, '(2x,f10.3)') Ez_drift_Vm
        WRITE (9, '("--dddddd.ddd- B_x (Gauss), not to be used in the equation of motion ---------")')
        WRITE (9, '(2x,f10.3)') Bx_drift_Gs

     END IF

  END IF

  CLOSE (9, STATUS = 'KEEP')
  
  piconst = 2.0_8 * asin(1.0_8)
  v_Te_ms     = SQRT(2.0_8 * T_e_eV * e_Cl / m_e_kg) !! defined with a factor of 2 !!
  W_plasma_s1 = SQRT(N_plasma_m3 * e_Cl**2 / (eps_0_Fm * m_e_kg))
!  W_cycl_x_s1 = e_Cl * B_x_ext_Gs * 0.0001_8 / m_e_kg
!  W_cycl_y_s1 = e_Cl * B_y_ext_Gs * 0.0001_8 / m_e_kg
  L_debye_m   = v_Te_ms / W_plasma_s1

  ! If the micron_flag is non-zero, then the cell size is specified directly.
  ! Be aware! If micron_flag is greater than zero than 
  !                       delta_x_m [ps] = N_max_vel*10^(n-1)
  !           If micron_flag is less than zero than
  !                       delta_x_m [ps] = N_max_vel*10^(n)
  if (micron_flag.ge.1) then
    delta_x_m = 1.e-6 * dble(N_of_cells_debye) * 10.**(micron_flag - 1)
    r_cells_debye =  L_debye_m / delta_x_m
  else if (micron_flag.lt.0) then
    delta_x_m = 1.e-6 * dble(N_of_cells_debye) * 10.**(micron_flag)
    r_cells_debye =  L_debye_m / delta_x_m
  else ! If micron_flag is zero, then the number of cells per Debye length is specified directly.
    r_cells_debye = dble(N_of_cells_debye)
    delta_x_m = L_debye_m / dble(N_of_cells_debye)
  end if

  ! If the picosec_flag is non-zero, then the time step is specified directly.
  ! Be aware! If picosec_flag is greater than zero than 
  !                       delta_t_s [ps] = N_max_vel*10^(n-1)
  !           If picosec_flag is less than zero than
  !                       delta_t_s [ps] = N_max_vel*10^(n)
  if (picosec_flag.ge.1) then
    delta_t_s = 1.e-12 * dble(N_max_vel) * 10.**(picosec_flag - 1) !to specify delta_t explicitly
    r_max_vel = delta_x_m / (v_Te_ms * delta_t_s)
  else if (picosec_flag.lt.0) then
    delta_t_s = 1.e-12 * dble(N_max_vel) * 10.**(picosec_flag) !to specify delta_t explicitly
    r_max_vel = delta_x_m / (v_Te_ms * delta_t_s)
  else
    delta_t_s = delta_x_m / (dble(N_max_vel) * v_Te_ms)
    r_max_vel = dble(N_max_vel)
  end if

  IF ((R_ext_ohm.GT.0.0_8).AND.(S_electrode_cm2.GT.0.0_8)) THEN
     factor_SR = (delta_x_m * delta_t_s) / (1.0d-4 * S_electrode_cm2 * R_ext_ohm * eps_0_Fm)
  ELSE
     factor_SR = 0.0_8
  END IF

!  factor_sigma = 1.0_8 / (e_Cl * N_in_macro) 
  factor_sigma = 1.0_8/(e_Cl * N_plasma_m3 * delta_x_m) * dble(N_of_particles_cell)
  factor_j     = delta_t_s * factor_sigma

  E_scl_Vm    = m_e_kg * v_Te_ms * W_plasma_s1 / e_Cl 
  U_scl_V     = E_scl_Vm * delta_x_m
  N_scl_m3    = N_plasma_m3 / DBLE(N_of_particles_cell)  
  V_scl_ms    = v_Te_ms / DBLE(N_box_vel)

  U_ext = U_ext_V / U_scl_V

  Max_T_cntr  = 1.0d-9 * T_sim_ns / delta_t_s                               ! Number of time steps

! if the plasma length is negative then it is given in beam wavelengths, assumed to be integer !!!
  IF (L_plasma_m.LT.0.0_8) THEN  
     IF (Rank_of_process.EQ.0) PRINT '(/2x,"The plasma length was given in units of electron beam wavelength")'
! to determine the length of the system as integer number of the wavelengths we will need the beam energy. Read it from the beam file.
     INQUIRE (FILE = 'ssc_ebeam.dat', EXIST = exists)
     IF(exists) THEN
        OPEN (9, FILE = 'ssc_ebeam.dat')
        READ (9, '(A77)') buf ! ********************** PERIODICITY OF PARTICLES MOTION **********************")')
        READ (9, '(A77)') buf ! ---           0 = turned off                                              ---")')
        READ (9, '(A77)') buf ! -------d----- 1 = turned on -------------------------------------------------")')
        READ (9, '(7x,i1)') dummy_i
        READ (9, '(A77)') buf ! ****************************** BEAM IN PLASMA *******************************")')
        READ (9, '(A77)') buf ! ---           0 = turned off                                              ---")')
        READ (9, '(A77)') buf ! ---           1 = monoenergetic beam                                      ---")')
        READ (9, '(A77)') buf ! -------d----- 2 = warm beam -------------------------------------------------")')
        READ (9, '(7x,i1)') dummy_i
        READ (9, '(A77)') buf ! --dddddd.ddd- Beam appears at this moment, [ns] (timesteps if < 0) ----------")')
        READ (9, '(2x,f10.3)') dummy_r
        READ (9, '(A77)') buf ! ============================= BEAM PARAMETERS ===============================")')
        READ (9, '(A77)') buf ! --dddddd.ddd- Energy, [eV] --------------------------------------------------")')
        READ (9, '(2x,f10.3)') Energy_ebeam_eV
        CLOSE (9, STATUS = 'keep')   !## this line was added on 2008-10-11
     ELSE
        IF  (Rank_of_process.EQ.0) THEN
           PRINT '(/2x,"ERROR: File with the name ssc_ebeam.dat not found ...")'
           PRINT  '(2x,"The system lengths was given in the wavelengthes, but the beam energy is unknown !")'
           PRINT  '(2x,"Program will be terminated now :(")'
        END IF
        STOP
     END IF
!     L_plasma_deb = - L_plasma_m
     L_plasma_deb = - L_plasma_m * 6.2831853_8 * SQRT(Energy_ebeam_eV / T_e_eV)  ! ####### ####### was like in the previous line ##### #####        
                                                ! ####### ####### now we assume that "-L_plasma_m" is the number of wavelengths ##### #####
  ELSE
     L_plasma_deb = L_plasma_m / L_debye_m
  END IF

!  N_cells  = L_plasma_deb * N_of_cells_debye      ! Number of cells is calculated    
  N_cells  = L_plasma_deb * r_cells_debye      ! Number of cells is calculated      ! ####### ####### was like in the previous line ##### #####

  N_nodes  = N_cells + 1                          ! Number of nodes is calculated

  CALL DETERMINE_NUMBER_OF_PARTICLES
!  N_of_part_initial = N_of_particles_cell * N_cells   

  L_plasma_m = N_cells * delta_x_m

! Produce filenames for checkpoints
  WRITE (procnumber_txt, '(i3)') Rank_of_process    ! Conversion of integer Rank_of_process to character*3 procnumber_txt
  procnumber_txt = ADJUSTL(TRIM(procnumber_txt))    ! Align to the left in order to
  proc_blnks = 3 - LEN_TRIM(procnumber_txt)         ! calculate the number of blanks in string number_txt;
  procnumber_txt = ADJUSTR(procnumber_txt)          ! then align to the right and fulfill
  procnumber_txt(1:proc_blnks) = '000'              ! substitution of leading blanks by '0'
  check_g_filename      = 'proc_PPP_g.check'
  check_g_filename(6:8) = procnumber_txt
  IF (Rank_of_process.GT.0) THEN
     check_e_filename      = 'proc_PPP_e.check'
     check_e_filename(6:8) = procnumber_txt
     check_i_filename      = 'proc_PPP_i.check'
     check_i_filename(6:8) = procnumber_txt
  END IF

  IF (M_i_amu.LT.0.0_8) THEN
     IF (Rank_of_process.EQ.0) THEN 
        PRINT '(/2x,"Ion mass is ",f10.3," electron masses or ",e12.5," kg")', Rank_of_process, ABS(M_i_amu), ABS(M_i_amu) * m_e_kg / amu_kg
     END IF
     M_i_amu = ABS(M_i_amu) * m_e_kg / amu_kg
  END IF

! set the drift velocity for the initial distribution
  IF ((Ez_drift_Vm.EQ.0.0_8).OR.(Bx_drift_Gs.EQ.0.0_8)) THEN
     Vy_e_drift = 0.0_8
  ELSE
     Vy_e_drift = N_box_vel * (10000.0_8 * Ez_drift_Vm / Bx_drift_Gs) / v_Te_ms
  END IF
  IF (Rank_of_process.EQ.0) THEN 
     PRINT '(/2x,"The E x B drift velocity used in electron distribution is ",e12.5," m/s")', Vy_e_drift * v_Te_ms / N_box_vel
  END IF

! REAL(8) KVx          ! Coefficients for equations of motion (X-PRE-FIN-moving)                         
  KVx = 1.0_8 / (r_max_vel * DBLE(N_box_vel))

! REAL(8) K_Q          ! Coefficient used in Poisson's equation
!  K_Q = 1.0_8 / (DBLE(N_of_cells_debye) * DBLE(N_of_particles_cell))
  K_Q = 1.0_8 / (r_cells_debye * DBLE(N_of_particles_cell))

!*** for dielectric layer at -d < x < 0:
  eps_param = eps_layer * delta_x_m / d_m

!  REAL(8) A_123_1      ! factor for coefficients before Ax^{n-1} for equations of motion (PRE-acceleration)
!!  A_123_1 =  0.5_8 * ( DBLE(N_box_vel) / (DBLE(N_max_vel) * DBLE(N_of_cells_debye)) )
  A_123_1 =  0.5_8 * ( DBLE(N_box_vel) / (r_max_vel * r_cells_debye) ) 

!  REAL(8) factor_KxEx  ! factor for coefficients before Ex^{n+1} for final X correction
  factor_KxEx = 1.0_8 / (r_max_vel * DBLE(N_box_vel))

  CALL PREPARE_CALCULATION_OF_MAG_FIELD_PROFILES

  CALL CONFIGURE_MESH_ARRAYS
  CALL SETZEROS_MESH_ARRAYS

  CALL SETVALUES_SPECIES_ARRAYS  ! moved this subroutine here since K_Xi is now allocated in CONFIGURE_MESH_ARRAYS

! initialize particles
  IF (Restore_from_checkpoint.EQ.0) THEN    ! for ordinary initialization

! initialize the seed values for random numbers generators
     IF (Rank_of_process.EQ.0) THEN
        call random_init
        call random_seed(put=seed)
        call sgrnd(seed(1))
!!        PRINT '(/2x,"Process ",i3," : Seed for the random number generator: ",i12)', Rank_of_process, I_random_seed
        PRINT '(/2x,"Process ",i3," : Seed for the random number generator: ",i12)', Rank_of_process, seed(1) 
        DO i = 1, N_of_processes - 1
           itmp = int(2.d0**31 * grnd())
           if (itmp.le.0) itmp = (2**31 - 1) + itmp  
           CALL MPI_SEND(itmp, 1, MPI_INTEGER, i, 101, MPI_COMM_WORLD, ierr)
        END DO
     ELSE
        CALL MPI_RECV(I_random_seed, 1, MPI_INTEGER, 0, 101, MPI_COMM_WORLD, stattus, ierr)
        call random_init
        seed(:) = I_random_seed
        call random_seed(put=seed)
        call sgrnd(seed(1))
        PRINT '(2x,"Process ",i3," : received the random numbers generator seed: ",i12)', Rank_of_process, I_random_seed
     END IF
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

! service data -------------------------------
     N_part_client_spec = 0

! determine numbers of particles, which will be assigned to the clients
     IF (N_of_processes.GT.1) THEN
        itmp = N_of_part_initial / (N_of_processes - 1)
        DO i = 1, N_of_processes - 1
           N_part_client_spec(i, 1) = itmp               ! for electrons
        END DO
        IF (itmp.LT.N_of_part_initial) N_part_client_spec(N_of_processes - 1, 1) = N_of_part_initial - itmp * (N_of_processes - 2)    ! if N_of_processes = 2
     END IF
     IF (N_spec.ge.2) N_part_client_spec(1:(N_of_processes-1), 2) = N_part_client_spec(1:(N_of_processes-1), 1)     ! for ions, if they are accounted
!     IF (N_spec.ge.3) N_part_client_spec(1:(N_of_processes-1), 3) = N_part_client_spec(1:(N_of_processes-1), 1)
! N_part(1:N_spec)     Numbers of particles
     N_part = 0
     IF (Rank_of_process.EQ.0) THEN                    ! for the server
        DO s = 1, N_spec
           N_part(s) = N_of_part_initial      
        END DO
     ELSE                                              ! for the client
        DO s = 1, N_spec
           N_part(s) = N_part_client_spec(Rank_of_process, s)
        END DO
     END IF

! Length(1:N_spec)     Lengthes of allocated arrays in blocks (they are not allocated yet!!!)
     DO s = 1, N_spec
        Length(s) = N_part(s) + INT(0.01 * REAL(N_part(s))) ! all processes
     END DO

     CALL CONFIGURE_PARTICLE_DYNAM_ARRAYS
     CALL INITIATE_PARTICLE_DYNAM_ARRAYS

     Q_left  = 0
     Q_right = 0
     Q_ext   = 0.0_8
     Start_T_cntr  = 0

     full_Q_left  = 0.0_8
     full_Q_right = 0.0_8
    
     prev_Q_left  = 0.0_8   !
     prev_Q_right = 0.0_8   !
     
  ELSE ! if we must read data from checkpoint files

     IF (Rank_of_process.EQ.0) THEN

        OPEN (9, FILE = check_g_filename)
        READ (9, '(2x,i14,2x,i14)') Start_T_cntr,     I_random_seed                   
        READ (9, '(2x,i14,2x,i14)') Start_diag_Tcntr, current_snap
        READ (9, '(2x,i14,2x,i14)') N_part(1),        N_part(2)
        READ (9, '(2x,i14,2x,i14)') Q_left,           Q_right
        READ (9, '(2x,e14.7,2x,e14.7,2x,e14.7)') full_Q_left, full_Q_right, Q_ext
        READ (9, '(2x,i14,2x,i14)')     N_of_saved_records, text_output_counter
        READ (9, '(2x,e14.7,2x,e14.7)') Energy_full_eV,     Init_energy_full_eV
        READ (9, '(2x,e14.7,2x,e14.7)') Energy_pot_eV,      Energy_heat_eV
        READ (9, '(2x,e14.7,2x,e14.7)') Energy_kin_eV(1),   Energy_kin_eV(2)
        READ (9, '(2x,e14.7,2x,e14.7)') Energy_wall_eV(1),  Energy_wall_eV(2)
        READ (9, '(2x,e14.7,2x,e14.7)') Energy_emit_eV(1),  Energy_emit_eV(2)
        READ (9, '(2x,e14.7,2x,e14.7)') Energy_coll_eV(1),  Energy_coll_eV(2)
        READ (9, '(2x,e14.7,2x,e14.7)') prev_Q_left,      prev_Q_right
!!        Q_ext = 0.0_8
       CLOSE (9, STATUS = 'KEEP')

       call random_init !*** "seed" allocated there
       seed(:) = I_random_seed
       call random_seed(put=seed)
       call sgrnd(seed(1))

     ELSE

        OPEN (9, FILE = check_g_filename)
        READ (9, '(2x,i14,2x,i14)') Start_T_cntr,     I_random_seed                   
        READ (9, '(2x,i14,2x,i14)') Start_diag_Tcntr, current_snap
        READ (9, '(2x,i14,2x,i14)') N_part(1),        N_part(2)
        READ (9, '(2x,i14,2x,i14)') Q_left,           Q_right
        CLOSE (9, STATUS = 'KEEP')

        call random_init
        seed(:) = I_random_seed
        call random_seed(put=seed)
        call sgrnd(seed(1))

        DO s = 1, N_spec
           Length(s) = N_part(s) + INT(0.01 * REAL(N_part(s)))
        END DO

        CALL CONFIGURE_PARTICLE_DYNAM_ARRAYS
        CALL READ_PARTICLE_DYNAM_ARRAYS

     END IF

     IF (Restore_from_checkpoint.EQ.2) Start_T_cntr = 0 

  END IF

  IF (Rank_of_process.EQ.0) THEN 
     PRINT '(/2x,"Electron Langmuir frequency                      : ",e10.3, " s^-1")', W_plasma_s1
     PRINT  '(2x,"     Ion Langmuir frequency                      : ",e10.3, " s^-1")', W_plasma_s1 * SQRT(m_e_kg / (M_i_amu * amu_kg))
!
!     PRINT  '(2x,"Electron cyclotron frequency                     : ",e10.3, " s^-1")', SQRT(W_cycl_x_s1**2 + W_cycl_y_s1**2)
!     PRINT  '(2x,"     Ion cyclotron frequency                     : ",e10.3, " s^-1")', SQRT(W_cycl_x_s1**2 + W_cycl_y_s1**2) * m_e_kg / (M_i_amu * amu_kg)
!
     PRINT  '(2x,"The electron thermal velocity                    : ",e10.3, " m/s")',  v_Te_ms
     PRINT  '(2x,"The electron Debye length                        : ",e10.3, " m")',    L_debye_m 

     PRINT '(/2x,"After mesh distribution the plasma length is     : ",f9.3," cm")', L_plasma_m * 100.0_8
!     PRINT  '(2x,"The new plasma length is                         : ",i6," (debye len.)")', (N_cells / N_of_cells_debye)
      PRINT  '(2x,"The new plasma length is                         : ",i6," (debye len.)")', int(N_cells / r_cells_debye)
     PRINT  '(2x,"Total number of plasma nodes is                  : ",i9)', N_nodes
     PRINT  '(2x,"The cell size is                                 : ",f9.3," mkm")', delta_x_m * 1.0e6
     IF (Restore_from_checkpoint.EQ.0) THEN
        PRINT  '(2x,"Initial total number of particles of one species : ",i9)', N_of_part_initial
     ELSE
        PRINT  '(2x,"Initial total number of particles (e)            : ",i9)', N_part(1)
        PRINT  '(2x,"Initial total number of particles (i)            : ",i9)', N_part(2)        
        PRINT  '(2x,"Charge assigned to the left wall                 : ",i9," charges of one particle")', Q_left
        PRINT  '(2x,"Charge assigned to the right wall                : ",i9," charges of one particle")', Q_right
     END IF
     PRINT  '(2x,"Number of time steps                             : ",i9)', Max_T_cntr
     PRINT  '(2x,"The time step is                                 : ",f9.3," ps")', delta_t_s * 1.0e12

     IF (rf_on) THEN
        PRINT  '(/2x,"The RF frequency is                              : ",e10.3," Hz")', f_rf_Hz
        PRINT  '(2x,"The RF amplitude is                              : ",e10.3," V")',  U_rf_V
        PRINT  '(2x,"The RF voltage will begin at                     : ",e10.3," s")', t_start_s
     END IF

     PRINT '(/2x,"Ratio of electron plasma period to the time step : ",f9.2)', 6.28318530718_8 / (W_plasma_s1 * delta_t_s)

     W_cycl_min_s1 = e_Cl * SQRT((Bx_gauss(0.0_8))**2 + (By_gauss(0.0_8))**2) * 0.0001_8 / m_e_kg
     W_cycl_max_s1 = W_cycl_min_s1
     i_Wc_min = 0
     i_Wc_max = 0

     DO i = 1, N_cells
        W_cycl_temp_s1 = e_Cl * SQRT((Bx_gauss(DBLE(i)))**2 + (By_gauss(DBLE(i)))**2) * 0.0001_8 / m_e_kg
        IF (W_cycl_temp_s1.LT.W_cycl_min_s1) THEN
           W_cycl_min_s1 = W_cycl_temp_s1
           i_Wc_min = i
        END IF
        IF (W_cycl_temp_s1.GT.W_cycl_max_s1) THEN
           W_cycl_max_s1 = W_cycl_temp_s1
           i_Wc_max = i
        END IF
     END DO

     PRINT '(/2x,"Minimal electron cyclotron frequency is ",e10.3," s^-1 at node ",i6," where x = ",f10.2," mm")', W_cycl_min_s1, i_Wc_min, i_Wc_min * delta_x_m * 1000.0_8
     IF (W_cycl_min_s1.NE.0.0_8) THEN
        PRINT  '(2x,"Here the ratio of electron gyro period to the time step is maximal and equals ",f9.2)', 6.28318530718_8 / (W_cycl_min_s1 * delta_t_s)
     END IF
     IF (N_spec.EQ.2) THEN
        PRINT '(2x,"   the minimal ion cyclotron frequency is ",e10.3," s^-1")', W_cycl_min_s1 / Ms(2)
        IF (W_cycl_min_s1.NE.0.0_8) THEN
           PRINT  '(2x,"   and the maximal ratio of ion gyro period to the time step is ",f11.2)', 6.28318530718_8 * Ms(2) / (W_cycl_min_s1 * delta_t_s)
        END IF
     END IF

     PRINT '(2x,"Maximal electron cyclotron frequency is ",e10.3," s^-1 at node ",i6," where x = ",f10.2," mm")', W_cycl_max_s1, i_Wc_max, i_Wc_max * delta_x_m * 1000.0_8
     IF (W_cycl_max_s1.NE.0.0_8) THEN
        PRINT  '(2x,"Here the ratio of electron gyro period to the time step is minimal and equals ",f9.2)', 6.28318530718_8 / (W_cycl_max_s1 * delta_t_s)
     END IF
     IF (N_spec.EQ.2) THEN
        PRINT '(2x,"   the maximal ion cyclotron frequency is ",e10.3," s^-1")', W_cycl_max_s1 / Ms(2)
        IF (W_cycl_max_s1.NE.0.0_8) THEN
           PRINT  '(2x,"   and the minimal ratio of ion gyro period to the time step is ",f11.2)', 6.28318530718_8 * Ms(2) / (W_cycl_max_s1 * delta_t_s)
        END IF
     END IF

  END IF

  call random_seed(get=seed)
END SUBROUTINE INITIATE_PARAMETERS 
 
!===================================================================================================
! Note, this subroutine is called in the INITIATE_PARAMETERS  
SUBROUTINE SETVALUES_SPECIES_ARRAYS 

  USE ParallelOperationValues
  USE CurrentProblemValues
  IMPLICIT NONE
  INTEGER s       ! species
  INTEGER i
  REAL(8) alfa_x, alfa_y, alfa_x2, alfa_y2, theta2, aa

! functions
  REAL(8) Bx_gauss
  REAL(8) By_gauss

! particle parameters -------------------------
  Ms  = 0.0_8
  Qs  = 0
  QMs = 0.0_8
  VT  = 0.0_8

  Ms(1)  =  1.0_8
  Qs(1)  = -1
  QMs(1) = -1.0_8
  VT(1)  =  dble(N_box_vel)
     
  DO s = 2, N_spec
     Ms(s)  = M_i_amu * amu_kg / m_e_kg
     Qs(s)  = 1          
     QMs(s) = dble(Qs(s)) / Ms(s)
     VT(s)  = SQRT( T_i_eV / (Ms(s) * T_e_eV) ) * dble(N_box_vel)
  END DO

  if (N_spec.eq.3) then !fast neutral species
     Ms(3) = Ms(2)
     Qs(3) = 0
     QMs(3) = 0.0_8
     VT(3) = VT(2)
  end if

! precalculated factors -----------------
  aa = r_max_vel * r_cells_debye

  DO s = 1, N_spec

     alfa_xy(s) = 0.5_8 * QMs(s) * (e_Cl * 0.0001_8 / m_e_kg) / (W_plasma_s1 * aa)

!!     A_123_3(s) = QMs(s) * (E_z_ext_Vm / E_scl_Vm) * (N_box_vel / aa)
     A_123_3(s) = 0.0_8 ! E_z_ext_Vm stores Ex in V/cm

     KvE_xyz(s) = 0.5_8 * QMs(s) * DBLE(N_box_vel) / aa

  END DO
  
  DO s = 1, N_spec
     DO i = 0, N_cells

        alfa_x = alfa_xy(s) * Bx_gauss(DBLE(i))
        alfa_y = alfa_xy(s) * By_gauss(DBLE(i))
        alfa_x2 = alfa_x**2
        alfa_y2 = alfa_y**2
        theta2 = alfa_x2 + alfa_y2
       
        K_Xi(i,s) = (0.25_8 * QMs(s) * Qs(s) / (DBLE(N_of_particles_cell) * aa**2)) * (1.0_8 + alfa_x2) / (1.0_8 + theta2) 

     END DO
  END DO

!  N_inject = 0  ! Number (counter) of injected particles of one species (1 = e, 2 = i) due to ionization or SEE

END SUBROUTINE SETVALUES_SPECIES_ARRAYS

!===================================================================================================
! Note, this subroutine is called in the INITIATE_PARAMETERS  
! Note, the new numbering assumes that node # 0       is the left plasma boundary  x = 0
!                                  and node # N_cells is the right plasma boundart x = L_plasma
SUBROUTINE CONFIGURE_MESH_ARRAYS 

  USE ParallelOperationValues
  USE CurrentProblemValues
  IMPLICIT NONE

  INTEGER ALLOC_ERR

!  PRINT '(4x,"Allocating mesh arrays ...")'

  ALLOCATE(EX(0:N_cells), STAT=ALLOC_ERR)
  
  ALLOCATE(GradEX(0:(N_cells-1)), STAT=ALLOC_ERR)
  
  ALLOCATE(F(0:N_cells), STAT=ALLOC_ERR)  
  
  ALLOCATE(Xi(0:N_cells), STAT=ALLOC_ERR)
  
  ALLOCATE(Q_stream(0:N_cells), STAT=ALLOC_ERR)
  
  ALLOCATE(Q_strm_spec(0:N_cells,1:N_spec), STAT=ALLOC_ERR)
  
  ALLOCATE(K_Xi(0:N_cells,1:N_spec), STAT=ALLOC_ERR)

!  PRINT '(4x,"Finished :)")'

END SUBROUTINE CONFIGURE_MESH_ARRAYS 

!===================================================================================================
!  Note, this subroutine is called in the INITIATE_PARAMETERS  
SUBROUTINE SETZEROS_MESH_ARRAYS

  USE CurrentProblemValues
  IMPLICIT NONE

!  PRINT '(4x,"Setting the values of mesh arrays equal to zero ...")'

  EX          = 0.0_8
  GradEX      = 0.0_8
  F           = 0.0_8
  Xi          = 0.0_8
  Q_stream    = 0.0_8
  Q_strm_spec = 0.0_8

! note, K_Xi is initialized elsewhere, in SETVALUES_SPECIES_ARRAYS

!  PRINT '(4x,"Finished :)")'

END SUBROUTINE SETZEROS_MESH_ARRAYS

!===================================================================================================
SUBROUTINE REMOVE_MESH_ARRAYS

  USE ParallelOperationValues
  USE CurrentProblemValues
  IMPLICIT NONE

  INTEGER DEALLOC_ERR

!  PRINT '("Deleting mesh arrays ...")'

  DEALLOCATE(EX, STAT=DEALLOC_ERR)
  DEALLOCATE(GradEX, STAT=DEALLOC_ERR)
  DEALLOCATE(F, STAT=DEALLOC_ERR)
  DEALLOCATE(Xi, STAT=DEALLOC_ERR)
  DEALLOCATE(Q_stream, STAT=DEALLOC_ERR)
  DEALLOCATE(Q_strm_spec, STAT=DEALLOC_ERR)
  DEALLOCATE(K_Xi, STAT=DEALLOC_ERR)

!  PRINT '(2x,"Finished :)")'

END SUBROUTINE REMOVE_MESH_ARRAYS

!===================================================================================================
SUBROUTINE CONFIGURE_PARTICLE_DYNAM_ARRAYS 

  USE ParallelOperationValues
  USE CurrentProblemValues
  IMPLICIT NONE
  INTEGER ALLOC_ERR
  INTEGER s 

!  PRINT '(4x,"Allocating particle dynamic arrays ...")'

  IF (.NOT.ALLOCATED(X_of_spec)) THEN
     ALLOCATE(X_of_spec(1:N_spec), STAT=ALLOC_ERR)
     IF(ALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in ALLOCATE X_of_spec !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF
  DO s = 1, N_spec
     NULLIFY(X_of_spec(s)%part)
     IF (.NOT.ASSOCIATED(X_of_spec(s)%part)) THEN
        ALLOCATE(X_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
        IF(ALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in ALLOCATE X_of_spec(s)%part !!!")', Rank_of_process
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
     END IF
  END DO
  
  IF (.NOT.ALLOCATED(VX_of_spec)) THEN
     ALLOCATE(VX_of_spec(1:N_spec), STAT=ALLOC_ERR)
     IF(ALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in ALLOCATE VX_of_spec !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF
  DO s = 1, N_spec
     NULLIFY(VX_of_spec(s)%part)
     IF (.NOT.ASSOCIATED(VX_of_spec(s)%part)) THEN
        ALLOCATE(VX_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
        IF(ALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in ALLOCATE VX_of_spec(s)%part !!!")', Rank_of_process
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
     END IF
  END DO

  IF (.NOT.ALLOCATED(VY_of_spec)) THEN
     ALLOCATE(VY_of_spec(1:N_spec), STAT=ALLOC_ERR)
     IF(ALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in ALLOCATE VY_of_spec !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF
  DO s = 1, N_spec
     NULLIFY(VY_of_spec(s)%part)
     IF (.NOT.ASSOCIATED(VY_of_spec(s)%part)) THEN
        ALLOCATE(VY_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
        IF(ALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in ALLOCATE VY_of_spec(s)%part !!!")', Rank_of_process
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
     END IF
  END DO

  IF (.NOT.ALLOCATED(VZ_of_spec)) THEN
     ALLOCATE(VZ_of_spec(1:N_spec), STAT=ALLOC_ERR)
     IF(ALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in ALLOCATE VZ_of_spec !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF
  DO s = 1, N_spec
     NULLIFY(VZ_of_spec(s)%part)
     IF (.NOT.ASSOCIATED(VZ_of_spec(s)%part)) THEN
        ALLOCATE(VZ_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
        IF(ALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in ALLOCATE VZ_of_spec(s)%part !!!")', Rank_of_process
           PRINT  '(2x,"Program will be terminated now :(")'
          STOP
        END IF
     END IF
  END DO

  IF (.NOT.ALLOCATED(AX_of_spec)) THEN
     ALLOCATE(AX_of_spec(1:N_spec), STAT=ALLOC_ERR)
     IF(ALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in ALLOCATE AX_of_spec !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF
  DO s = 1, N_spec
     NULLIFY(AX_of_spec(s)%part)
     IF (.NOT.ASSOCIATED(AX_of_spec(s)%part)) THEN
        ALLOCATE(AX_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
        IF(ALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in ALLOCATE AX_of_spec(s)%part !!!")', Rank_of_process
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
     END IF
  END DO

  IF (.NOT.ALLOCATED(Tag_of_spec)) THEN
     ALLOCATE(Tag_of_spec(1:N_spec), STAT=ALLOC_ERR)
     IF(ALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in ALLOCATE Tag_of_spec !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF
  DO s = 1, N_spec
     NULLIFY(Tag_of_spec(s)%part)
     IF (.NOT.ASSOCIATED(Tag_of_spec(s)%part)) THEN
        ALLOCATE(Tag_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
        IF(ALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in ALLOCATE Tag_of_spec(s)%part !!!")', Rank_of_process
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
     END IF
  END DO
  
  IF (.NOT.ALLOCATED(prev_VX_of_spec)) THEN
     ALLOCATE(prev_VX_of_spec(1:N_spec), STAT=ALLOC_ERR)
     IF(ALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in ALLOCATE prev_VX_of_spec !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF
  DO s = 1, N_spec
     NULLIFY(prev_VX_of_spec(s)%part)
     IF (.NOT.ASSOCIATED(prev_VX_of_spec(s)%part)) THEN
        ALLOCATE(prev_VX_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
        IF(ALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in ALLOCATE prev_VX_of_spec(s)%part !!!")', Rank_of_process
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
     END IF
  END DO

!  PRINT '(2x,"Finished :)")'

END SUBROUTINE CONFIGURE_PARTICLE_DYNAM_ARRAYS 

!===================================================================================================
SUBROUTINE INITIATE_PARTICLE_DYNAM_ARRAYS 

  use mpi
  USE ParallelOperationValues
  USE CurrentProblemValues
  IMPLICIT NONE

!  INCLUDE 'mpif.h'

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
!  LOGICAL flag

  INTEGER i, s, k_start, k_fin !, k

  integer dimofbufer
  real(8), allocatable :: bufer(:) !(1:1000000)
  integer alloc_err

  real(8) factor_x, factor_z

  PRINT '(4x,"Process ",i3," : Initiating particle dynamic arrays ...", i8,2x,i8)', Rank_of_process, N_part(1), N_part(2)

  IF (Rank_of_process.EQ.0) THEN  ! server

     dimofbufer = 1000
     do s = 1, N_spec
        do i = 1, N_of_processes - 1
           dimofbufer = MAX(dimofbufer, N_part_client_spec(i, s)) 
        end do
     end do
     allocate(bufer(1:dimofbufer), stat = alloc_err)

     CALL INITIAL_DISTRIB_IN_SPACE
     CALL INITIAL_DISTRIB_IN_VELOCITY

     DO i = 1, 3
        CALL SHUFFLE_ARRAYS
     END DO

! introduce anisotropy and drift for the electrons, having sampled from a spherically-symm. distribution
     factor_x = SQRT(Tx_e_eV / T_e_eV)
     factor_z = SQRT(Tz_e_eV / T_e_eV)
     DO i = 1, N_part(1)
        VX_of_spec(1)%part(i) = VX_of_spec(1)%part(i) * factor_x
        VY_of_spec(1)%part(i) = VY_of_spec(1)%part(i) * factor_z + Vy_e_drift
        VZ_of_spec(1)%part(i) = VZ_of_spec(1)%part(i) * factor_z
     END DO

     DO s = 1, N_spec

        k_start = 1
        k_fin   = 0
        DO i = 1, N_of_processes - 1

           k_fin   = k_fin + N_part_client_spec(i, s)

           bufer(1:N_part_client_spec(i, s)) = X_of_spec(s)%part(k_start:k_fin)
           CALL MPI_SEND(bufer, N_part_client_spec(i, s), MPI_DOUBLE_PRECISION, i, 1, MPI_COMM_WORLD, ierr) 

           bufer(1:N_part_client_spec(i, s)) = VX_of_spec(s)%part(k_start:k_fin)
           CALL MPI_SEND(bufer, N_part_client_spec(i, s), MPI_DOUBLE_PRECISION, i, 2, MPI_COMM_WORLD, ierr) 

           bufer(1:N_part_client_spec(i, s)) = VY_of_spec(s)%part(k_start:k_fin)
           CALL MPI_SEND(bufer, N_part_client_spec(i, s), MPI_DOUBLE_PRECISION, i, 3, MPI_COMM_WORLD, ierr) 

           bufer(1:N_part_client_spec(i, s)) = VZ_of_spec(s)%part(k_start:k_fin)
           CALL MPI_SEND(bufer, N_part_client_spec(i, s), MPI_DOUBLE_PRECISION, i, 4, MPI_COMM_WORLD, ierr) 

           k_start = k_fin + 1

        END DO

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

     END DO

     CALL REMOVE_PARTICLE_DYNAM_ARRAYS          ! server process does not need these data any more

  ELSE                                                        ! client process >>>>

     dimofbufer = 1000
     do s = 1, N_spec
        dimofbufer = MAX(dimofbufer, N_part(s)) 
     end do
     allocate(bufer(1:dimofbufer), stat = alloc_err)

     DO s = 1, N_spec

        CALL MPI_PROBE(0, 1, MPI_COMM_WORLD, stattus, ierr)
        
        CALL MPI_RECV(bufer, N_part(s), MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, stattus, ierr)
        X_of_spec(s)%part(1:N_part(s))   = bufer(1:N_part(s))
        
        CALL MPI_RECV(bufer, N_part(s), MPI_DOUBLE_PRECISION, 0, 2, MPI_COMM_WORLD, stattus, ierr)
        VX_of_spec(s)%part(1:N_part(s))  = bufer(1:N_part(s))
        
        CALL MPI_RECV(bufer, N_part(s), MPI_DOUBLE_PRECISION, 0, 3, MPI_COMM_WORLD, stattus, ierr)
        VY_of_spec(s)%part(1:N_part(s))  = bufer(1:N_part(s))
        
        CALL MPI_RECV(bufer, N_part(s), MPI_DOUBLE_PRECISION, 0, 4, MPI_COMM_WORLD, stattus, ierr)
        VZ_of_spec(s)%part(1:N_part(s))  = bufer(1:N_part(s))

        Tag_of_spec(s)%part(1:N_part(s)) = 0
        AX_of_spec(s)%part(1:N_part(s))  = 0.0_8

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

     END DO

  END IF

  deallocate(bufer, stat = alloc_err)

  PRINT '(2x,"Finished :)")'

END SUBROUTINE INITIATE_PARTICLE_DYNAM_ARRAYS 

!===================================================================================================
!
SUBROUTINE READ_PARTICLE_DYNAM_ARRAYS

  USE ParallelOperationValues
  USE CurrentProblemValues

  IMPLICIT NONE
  INTEGER n

  IF (Rank_of_process.EQ.0) RETURN

  OPEN (9, FILE = check_e_filename)
  PRINT '(4x,"Process ",i3," : Reading ",i7," particles (e) from file ",A15)', Rank_of_process, N_part(1), check_e_filename
  DO n = 1, N_part(1)
     READ (9, '(5(1x,e16.9),1x,i4)')  X_of_spec(1)%part(n), &
                                   & VX_of_spec(1)%part(n), & 
                                   & VY_of_spec(1)%part(n), &
                                   & VZ_of_spec(1)%part(n), &
                                   & AX_of_spec(1)%part(n), &
                                   & Tag_of_spec(1)%part(n)
  END DO
  CLOSE (9, STATUS = 'KEEP')
     
  IF (N_spec.EQ.2) THEN
     OPEN (9, FILE = check_i_filename)
     PRINT '(4x,"Process ",i3," : Reading ",i7," particles (i) from file ",A15)', Rank_of_process, N_part(2), check_i_filename
     DO n = 1, N_part(2)
        READ (9, '(5(1x,e16.9),1x,i4)')  X_of_spec(2)%part(n), &
                                      & VX_of_spec(2)%part(n), & 
                                      & VY_of_spec(2)%part(n), &
                                      & VZ_of_spec(2)%part(n), &
                                      & AX_of_spec(2)%part(n), &
                                      & Tag_of_spec(2)%part(n)
     END DO
     CLOSE (9, STATUS = 'KEEP')
  END IF

END SUBROUTINE READ_PARTICLE_DYNAM_ARRAYS

!===================================================================================================
SUBROUTINE REMOVE_PARTICLE_DYNAM_ARRAYS

  USE ParallelOperationValues
  USE CurrentProblemValues
  IMPLICIT NONE
  INTEGER DEALLOC_ERR

!  PRINT '("Deleting particle dynamic arrays ...")'

  IF (ALLOCATED(X_of_spec)) THEN
     DEALLOCATE(X_of_spec, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE X_of_spec !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  IF (ALLOCATED(VX_of_spec)) THEN
     DEALLOCATE(VX_of_spec, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE VX_of_spec !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  IF (ALLOCATED(VY_of_spec)) THEN
     DEALLOCATE(VY_of_spec, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE VY_of_spec !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  IF (ALLOCATED(VZ_of_spec)) THEN
     DEALLOCATE(VZ_of_spec, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE VZ_of_spec !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  IF (ALLOCATED(AX_of_spec)) THEN
     DEALLOCATE(AX_of_spec, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE AX_of_spec !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  IF (ALLOCATED(Tag_of_spec)) THEN
     DEALLOCATE(Tag_of_spec, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE Tag_of_spec !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  IF (ALLOCATED(prev_VX_of_spec)) THEN
     DEALLOCATE(prev_VX_of_spec, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE prev_VX_of_spec !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

!  PRINT '(2x,"Finished :)")'

END SUBROUTINE REMOVE_PARTICLE_DYNAM_ARRAYS

!=============================================================================================
! preallocating the linked list for electrons
!
SUBROUTINE PREALLOCATE_ELECTRON_LIST

  USE ParallelOperationValues
  USE CurrentProblemValues
  IMPLICIT NONE
  INTEGER ALLOC_ERR

  NULLIFY(Inj_electron)

  IF (.NOT.ASSOCIATED(Inj_electron)) THEN
     ALLOCATE(Inj_electron, STAT=ALLOC_ERR)
     IF (ALLOC_ERR.NE.0) THEN
        PRINT '(/2x,"Process ",i3," : Error in ALLOCATE Inj_electron !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  Current_electron => Inj_electron
  NULLIFY(Current_electron%next)

END SUBROUTINE PREALLOCATE_ELECTRON_LIST

!---------------------------------------------------------------------------------------------
! clear the linked list for electrons
!
SUBROUTINE CLEAR_ELECTRON_LIST

  USE ParallelOperationValues
  USE CurrentProblemValues
  IMPLICIT NONE
  INTEGER DEALLOC_ERR

  TYPE (injected_particle), POINTER :: current, temp

  current => Inj_electron

  DO WHILE(ASSOCIATED(current))
     temp => current
     current => current%next
     DEALLOCATE(temp, STAT=DEALLOC_ERR)
     IF (DEALLOC_ERR.NE.0) THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE Inj_electron !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END DO

END SUBROUTINE CLEAR_ELECTRON_LIST

!=============================================================================================
! preallocating the linked list for ions
!
SUBROUTINE PREALLOCATE_ION_LIST

  USE ParallelOperationValues
  USE CurrentProblemValues
  IMPLICIT NONE
  INTEGER ALLOC_ERR

  NULLIFY(Inj_ion)

  IF (.NOT.ASSOCIATED(Inj_ion)) THEN
     ALLOCATE(Inj_ion, STAT=ALLOC_ERR)
     IF (ALLOC_ERR.NE.0) THEN
        PRINT '(/2x,"Process ",i3," : Error in ALLOCATE Inj_ion !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  Current_ion => Inj_ion
  NULLIFY(Current_ion%next)

END SUBROUTINE PREALLOCATE_ION_LIST

!---------------------------------------------------------------------------------------------
! clear the linked list for ions
!
SUBROUTINE CLEAR_ION_LIST

  USE ParallelOperationValues
  USE CurrentProblemValues
  IMPLICIT NONE
  INTEGER DEALLOC_ERR

  TYPE (injected_particle), POINTER :: current, temp

  current => Inj_ion

  DO WHILE(ASSOCIATED(current))
     temp => current
     current => current%next
     DEALLOCATE(temp, STAT=DEALLOC_ERR)
     IF (DEALLOC_ERR.NE.0) THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE Inj_ion !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END DO

END SUBROUTINE CLEAR_ION_LIST

!=============================================================================================
! This subroutine the number of particles of a single species, which is necessary to 
! reproduce the desired density profile.
! 
SUBROUTINE DETERMINE_NUMBER_OF_PARTICLES

  USE ParallelOperationValues
  USE CurrentProblemValues
  IMPLICIT NONE
  REAL(8) density
  REAL(8) x_left
  REAL(8) x_right
  REAL(8) dx

  INTEGER n_int
  INTEGER i

  REAL(8) integral
  REAL(8) x

  N_of_part_initial = N_of_particles_cell * N_cells
  
  n_int   = N_of_part_initial * 20
  x_left  = 0.5_8 / N_of_particles_cell !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! step aside from the walls !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  x_right = N_cells - x_left
  dx = (x_right - x_left) / n_int     

! Calculate the integral of the function "density(z)", which describes the density profile.
  integral = 0.0_8
  x = x_left
  DO i = 1, n_int 
     x = x_left + i * dx            
     integral = integral + density(x - 0.5_8 * dx)
  END DO

  N_of_part_initial = N_of_particles_cell * (integral * dx)

END SUBROUTINE DETERMINE_NUMBER_OF_PARTICLES


!=============================================================================================
! This subroutine distributed the electrons and ions in space according to the density profile.
! 
SUBROUTINE INITIAL_DISTRIB_IN_SPACE

  USE ParallelOperationValues
  USE CurrentProblemValues
  IMPLICIT NONE
  REAL(8) density
  REAL(8) x_left
  REAL(8) x_right
  REAL(8) dx
  REAL(8) totintegr
  INTEGER n_int
  INTEGER i, j, n
!  INTEGER indx

  REAL(8), ALLOCATABLE :: integral(:)
  REAL(8), ALLOCATABLE :: x(:)
  INTEGER ALLOC_ERR, DEALLOC_ERR
  INTEGER i_prev     ! These integer variables are used to trace the following situation: 
  INTEGER i_new      ! a number ABC{D}.**** increases to ABC{D+1}.**** (ABCD - decimal digits) 

  PRINT '(6x,"Process ",i3," : Performing initial distribution of particles in space ...")', Rank_of_process 

! for testing purposes
! n= 1
!   do i = 1, N_part(n)
!      X_of_spec(n)%part(i) = 0.1_8 * REAL(N_cells)
!   end do
!end do
!do n= 2, N_spec
!   do i = 1, N_part(n)
!      X_of_spec(n)%part(i) = 0.9_8 * REAL(N_cells)
!   end do
!end do
!return
  
! Note, first particle will be always distributed to x_left, therefore integral is used to distribute (N_of_part_initial - 1) particle
  n_int   = (N_of_part_initial - 1)  * 10 ! was 20
  x_left  = 0.5_8 / N_of_particles_cell !+ 20.0_8 !0.0_8  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! step aside from the walls !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  x_right = N_cells - x_left
  dx = (x_right - x_left) / n_int     
! Allocate array of values of the integral of density profile function
  IF (.NOT.ALLOCATED(integral)) THEN
     ALLOCATE(integral(0:n_int), STAT=ALLOC_ERR)
     IF(ALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in ALLOCATE density integral by x !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF
! Allocate array of corresponding values of X-coordinate
  IF (.NOT.ALLOCATED(x)) THEN
     ALLOCATE(x(0:n_int), STAT=ALLOC_ERR)
     IF(ALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in ALLOCATE x(0:",i9,") !!!")', Rank_of_process, n_int
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF
! Calculate and store the values of integral of the function "density(z)",
! which describes the density profile.
  integral(0) = 0.0_8
  x(0) = x_left
  DO i = 1, n_int 
     x(i) = x_left + i * dx            
     integral(i) = integral(i-1) + density(x(i) - 0.5_8 * dx)
  END DO
! Normalize the values of the integral
  totintegr = integral(n_int)
  integral = dble(N_of_part_initial) * integral / totintegr
! Distribute the particles
  n = 1                       ! particles counter
  i_prev = -1                 ! set so to distribute the very first particle
  DO i = 0, n_int 
     i_new = int(integral(i))
! If the integral value increases by 1 we place the particles (electron and ion) at this point.
     IF ((i_new - i_prev).GE.1) THEN
! Here we assign the same coordinates to the current electrons (j=1) and ions (j=2)
        DO j = 1, N_spec   
           X_of_spec(j)%part(n) = x(i) 
        END DO
        i_prev = n
        n = n + 1      ! increase the particles counter
        IF ((n.GT.N_of_part_initial).AND.(i.LT.n_int)) THEN
           PRINT '(2x,"Process ",i3," : Steps left (to finish the space integration): ",i9)', Rank_of_process, n_int - i
           EXIT
        END IF
     END IF
!!     i_prev = integral(i)
  END DO
  PRINT '(/4x,"Process ",i3," : ",i8," particles of each species were distributed in space")', Rank_of_process, n-1
! Check the possible situation when some particles would be not included into the distribution.
  IF (n.LE.N_of_part_initial) THEN
     PRINT '(/2x,"Process ",i3," : Error! The last particles were not distributed in space: ",i9)', &
                                                                        &Rank_of_process, N_of_part_initial - n + 1
     PRINT  '(2x,"Program will be terminated now :(")'
     STOP
  END IF
! Deallocate arrays 
  IF (ALLOCATED(integral)) THEN
     DEALLOCATE(integral, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE density integral over x !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF
  IF (ALLOCATED(x)) THEN
     DEALLOCATE(x, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE x !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  PRINT '(2x,"Finished :)")'

END SUBROUTINE INITIAL_DISTRIB_IN_SPACE

!==========================================================
REAL(8) FUNCTION density(x)

  USE CurrentProblemValues, ONLY : N_cells, N_distrib_m3, N_plasma_m3, Density_flag

  IMPLICIT NONE
  
  REAL(8) x
! We assume the constant profile of density distribution

!  density = 1.0_8

  IF (N_distrib_m3.EQ.N_plasma_m3) THEN
     density = 1.0_8
  ELSE   
     density = N_distrib_m3 / N_plasma_m3
!     IF (x.GT.(0.95_8*N_cells)) density = 0.0_8
!     IF (x.LT.(0.05_8*N_cells)) density = 0.0_8
   END IF
 
!  density = 0.01_8

!  density = 0.25_8

  IF (Density_flag.EQ.1) THEN
     density = density * (-6.0_8 * (x / N_cells - 0.5_8)**2 + 2.0_8) * (2.0_8 / 3.0_8) ! 2/3 normalization factor
  END IF

!density  = 0.0_8  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!if ((x.ge.0.0_8).and.(x.le.100.0_8)) density = 1.0_8    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!if((x.lt.100.0_8).or.(x.gt.505.0_8)) density = 0.1_8
!if((x.gt.600.0_8)) density = 0.1_8

END FUNCTION density

!==========================================================
SUBROUTINE INITIAL_DISTRIB_IN_VELOCITY

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE mt19937
  IMPLICIT NONE
  REAL(8) maxwell
  REAL(8) v_min
  REAL(8) v_max
  REAL(8) dv, y, gauss, s0, s2  
  REAL(8) totintegr
  INTEGER n_int
  INTEGER i, j
  INTEGER n           ! the counter of the particles
!  INTEGER indx
  real r1, r2, costheta, sintheta, phi
  REAL(8), ALLOCATABLE :: integral(:)
  REAL(8), ALLOCATABLE :: v(:)
  INTEGER ALLOC_ERR, DEALLOC_ERR
  INTEGER i_prev     ! These integer variables are used to trace the following situation: 
  INTEGER i_new      ! a number ABC{D}.**** increases to ABC{D+1}.**** (ABCD - decimal digits) 
  
  PRINT '(6x,"Process ",i3," : Performing initial distribution of particles over velocity ...")', Rank_of_process

  call random_seed(put=seed)

!! for testing purposes only (ions do not matter here, we will watch electrons only) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!do j = 1, N_spec
!   do n = 1, N_part(j)
!      VX_of_spec(j)%part(n) = 0.0_8 !VT(j) * sqrt(initial_electron_energy_eV / T_e_eV) * COS(initial_angle_of_incidence_deg * 3.141592654_8 / 180.0_8)
!!if (mod(n,2).EQ.0) VX_of_spec(j)%part(n) = -VX_of_spec(j)%part(n)   ! we want two counterpropagating beams
!      VY_of_spec(j)%part(n) = 0.0_8 !VT(j) * sqrt(initial_electron_energy_eV / T_e_eV) * SIN(initial_angle_of_incidence_deg * 3.141592654_8 / 180.0_8)
!      VZ_of_spec(j)%part(n) = 0.0_8
!   end do
!end do 
!return
  
  n_int = (N_of_part_initial-1) * 20  !10
! Allocate array of values of the integral of density profile function
  IF (.NOT.ALLOCATED(integral)) THEN
     ALLOCATE(integral(0:n_int), STAT=ALLOC_ERR)
     IF(ALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in ALLOCATE integral over v !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF
! Allocate array of corresponding values of velocity
  IF (.NOT.ALLOCATED(v)) THEN
     ALLOCATE(v(0:n_int), STAT=ALLOC_ERR)
     IF(ALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in ALLOCATE v !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF
! For each type (j) of particles
  DO j = 1, N_spec
!!     v_min = -VT(j) * dble(N_max_vel_distrib) 
     v_min = 0.
     v_max =  VT(j) * dble(N_max_vel_distrib) 
     dv = (v_max - v_min) / dble(n_int)
! Calculate the initial value of integral of the function maxwell(z),
! which describes the velocity distribution density profile.
! It is necessary for normalization of the density function.
     integral(0) = 0.0_8
     v(0) = -0.5_8 * dv
     DO i = 1, n_int 
        v(i) = v_min + (dble(i)-0.5_8) * dv
        y = v(i) / VT(j)
        gauss = dexp(-y**2)
        s0 = 0.5 * sqrt(piconst) * erf(y)
        s2 = 0.5 * (s0 - y * gauss)
        integral(i) = s2
!        integral(i) = integral(i-1) + maxwell(v(i), VT(j)) * (v(i)**2)
     END DO
! Normalize the values of the integral
     totintegr = integral(n_int)
     integral = dble(N_of_part_initial + 1) * integral / totintegr
! Distribute the particles
     n = 1 !** index of the next particle to be loaded
     i_prev = 0
     DO i = 1, n_int 
        i_new = int(integral(i))
! If the integral value increases by 1 we set the particle velocity to the current value of v.
         IF ((i_new - i_prev).GE.1) THEN   
!!           call random_number(r1)
!!           call random_number(r2)
           r1 = grnd()
           costheta = 2. * r1 - 1.
           sintheta = 2. * sqrt(r1 * (1. - r1))
           phi = 2. * piconst * grnd()
           VX_of_spec(j)%part(n) = v(i) * costheta
           VY_of_spec(j)%part(n) = v(i) * sintheta * cos(phi)
           VZ_of_spec(j)%part(n) = v(i) * sintheta * sin(phi)
           i_prev = n
           n = n + 1            ! increase the particles counter
           IF ((n.GT.N_of_part_initial).AND.(i.LT.n_int)) THEN
              PRINT '(2x,"Process ",i3," : Steps left (to finish the velocity integration): ",i9)', Rank_of_process, n_int - i
              EXIT
           END IF
        END IF
     END DO
     PRINT '(4x,"Process ",i3," : ",i8," particles of species ",i1," were distributed over velocity")',Rank_of_process, n-1, j
! Check the possible situation when some particles would be not included into the distribution.
     IF (n.LE.N_of_part_initial) THEN
        PRINT '(/2x,"Process ",i3," : Warning! The last particles were not distributed over velocity: ",i9)', &
                                                                & Rank_of_process, N_of_part_initial - n + 1
        DO i = n, N_of_part_initial
           VX_of_spec(j)%part(i) = 0.0_8
           VY_of_spec(j)%part(i) = 0.0_8
           VZ_of_spec(j)%part(i) = 0.0_8
        END DO
        PRINT '(2x,"These particles were assigned zero velocities before shuffling")'
!        PRINT  '(2x,"Program will be terminated now :(")'
!        STOP
     END IF
  END DO
! Deallocate arrays 
  IF (ALLOCATED(integral)) THEN
     DEALLOCATE(integral, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE integral by v !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF
  IF (ALLOCATED(v)) THEN
     DEALLOCATE(v, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE v !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  PRINT '(2x,"Finished :)")'
  call random_seed(get=seed)
END SUBROUTINE INITIAL_DISTRIB_IN_VELOCITY

!==========================================================
REAL(8) FUNCTION maxwell(v, vt)

  IMPLICIT NONE
  REAL(8) v, vt
! We assume the isothermic maxwell plasma

!maxwell = 1.0_8
!if (abs(v).gt.vt) maxwell = 0.0_8

  maxwell = DEXP(-v*v/(vt*vt))   ! 1.0_8                          !!!!!!!!&^*^%&$%&^%R(&^)*(&*(&)(^

!  maxwell = DEXP(-4.0_8 * v*v/(vt*vt))   ! 1.0_8         ! ####### ####### we use the 4 times smaller temperature for distribution

!  maxwell = DEXP(-2.0_8 * v*v/(vt*vt))   ! 1.0_8         ! ####### ####### we use the 2 times smaller temperature for distribution

! if (abs(v).gt.vt) maxwell = 0.0_8


END FUNCTION maxwell

!==========================================================
! This subroutine MUST be called AFTER the subroutines
!  INITIAL_DISTRIBUTION_IN_SPACE
!  INITIAL_DISTRIBUTION_IN_VELOCITY
! in order to avoid the correlation
SUBROUTINE SHUFFLE_ARRAYS
  
  USE ParallelOperationValues
  USE CurrentProblemValues
  USE mt19937
  IMPLICIT NONE

  REAL RAN, R
  REAL(8) tmp_value
!  INTEGER itmp_value
  INTEGER i, s, random_i
!  INTEGER indx

  call random_seed(put=seed)

!!!return     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  PRINT '(6x,"Process ",i3," : Shuffling the particle dynamics arrays ...")', Rank_of_process

! First of all we shuffle the coordinates
  DO s = 1, N_spec
     DO i = 1, N_part(s) ! N_of_part_initial !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        tmp_value  =       X_of_spec(s)%part(i) 
! We use the unconditional DO loop here in order to avoid the situation when
! random_i < 1 or random_i > N_part
        DO
! The new random value will be retrieved until it will satisfy the desired range
! 1 <= random_i <= N_of_part_initial
           R= grnd()
           random_i = 1 + int(R * float(N_of_part_initial))
           IF ((random_i.GE.1).AND.(random_i.LE.N_of_part_initial)) THEN
              X_of_spec(s)%part(i)        = X_of_spec(s)%part(random_i)
              X_of_spec(s)%part(random_i) = tmp_value
              EXIT
           END IF
        END DO
     END DO
  END DO
! Shuffle the X-velocity
  DO s = 1, N_spec
     DO i = 1, N_part(s) ! N_of_part_initial !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        tmp_value = VX_of_spec(s)%part(i) 
        DO
           R= grnd()
           random_i = 1 + int(R * float(N_of_part_initial))
           IF ((random_i.GE.1).AND.(random_i.LE.N_of_part_initial)) THEN
              VX_of_spec(s)%part(i)        = VX_of_spec(s)%part(random_i)
              VX_of_spec(s)%part(random_i) = tmp_value
              EXIT
           END IF
        END DO
     END DO
  END DO
! Shuffle the Y-velocity
  DO s = 1, N_spec
     DO i = 1, N_part(s) ! N_of_part_initial !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        tmp_value = VY_of_spec(s)%part(i) 
        DO
           R= grnd()
           random_i = 1 + int(R * float(N_of_part_initial))
           IF ((random_i.GE.1).AND.(random_i.LE.N_of_part_initial)) THEN
              VY_of_spec(s)%part(i)        = VY_of_spec(s)%part(random_i)
              VY_of_spec(s)%part(random_i) = tmp_value
              EXIT
           END IF
        END DO
     END DO
  END DO
! Shuffle the Z-velocity
  DO s = 1, N_spec
     DO i = 1, N_part(s) ! N_of_part_initial !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        tmp_value  = VZ_of_spec(s)%part(i) 
        DO
           R= grnd()
           random_i = 1 + int(R * float(N_of_part_initial))
           IF ((random_i.GE.1).AND.(random_i.LE.N_of_part_initial)) THEN
              VZ_of_spec(s)%part(i)        = VZ_of_spec(s)%part(random_i)
              VZ_of_spec(s)%part(random_i) = tmp_value
              EXIT
           END IF
        END DO
     END DO
  END DO

!  PRINT '(2x,"Finished :)")'

  call random_seed(get=seed)
END SUBROUTINE SHUFFLE_ARRAYS













