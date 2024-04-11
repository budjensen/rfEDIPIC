!========================================================================================================================
!
module ParallelOperationValues

   integer Rank_of_process
   integer N_of_processes

end module ParallelOperationValues

!========================================================================================================================
!
module CurrentProblemValues

   integer, parameter :: N_spec = 2                   !! The number of species of particles used in simulations

   real(8), parameter :: e_Cl     = 1.602189d-19      !! Charge of single electron [Cl]
   real(8), parameter :: m_e_kg   = 9.109534d-31      !! Mass of single electron [kg]
   real(8), parameter :: eps_0_Fm = 8.854188d-12      !! The dielectric constant [F/m]
   real(8), parameter :: mu_0_Hm  = 1.256637d-6       !! The magnetic constant [H/m]
   real(8), parameter :: amu_kg = 1.660565d-27        !! atomic mass unit [kg]
!** properties of the dielectric layer positioned at -d < x < 0:
   real(8), parameter :: d_m = 0.01                   !! thickness of dielectric layer on top of the left electrode
   real(8), parameter :: eps_layer = 3.9_8
!!  real(8), parameter :: eps_layer = 1.e10 !*** to test implementation of the boundary condition
   real(8) eps_param

   real(8) R_ext_ohm       !! Resistance of resistor in external circuit [ohm]
   real(8) S_electrode_cm2 !! Area of an electrode [cm^2]
   real(8) E_z_ext_Vm      !! The external accelerating Z-electric field [V/m]
!  real(8) B_x_ext_Gs      ! The external radial X-magnetic field [Gauss]
!  real(8) B_y_ext_Gs      ! The external radial Y-magnetic field [Gauss]
   real(8) U_ext_V         !! The externally maintained potential difference between the metal end plates [V]
   real(8) f_rf_Hz         !! The frequency of the rf power source on the left wall (Hz)
   real(8) U_rf_V          !! The left wall rf power source voltage amplitude (V)
   real(8) t_start_s       !! The ignition time of the left wall's rf power source (s)

   logical rf_on           !! flag, turns on/off the rf power source on the left wall

   real(8) L_plasma_m      !! Plasma's length [m], aproximate, may be modified after mesh distribution
   real(8) N_plasma_m3     !! Plasma's density [m^-3], affects scaling
   real(8) M_i_amu         !! Mass of single ion [a.m.u.]
   real(8) T_e_eV          !! The electron temperature [eV], affects scaling
   real(8) T_i_eV          !! The ion temperature [eV]
   real(8) T_sim_ns        !! Duration of simulations [ns]
   real(8) piconst

   real(8) Tx_e_eV         !! x-electron temperature for anisotropic electron distribution [eV]
   real(8) Tz_e_eV         !! z-electron temperature for anisotropic electron distribution [eV]
   real(8) N_distrib_m3    !! Plasma density for electron distribution [m^-3], affects the number of particles, not the scaling
   real(8) Ez_drift_Vm     !! Electric field [V/m], for calculation of the ExB-drift y-velocity, not used in the equations of motion
   real(8) Bx_drift_Gs     !! Magnetic field [gauss], for calculation of the ExB-drift y-velocity, not used in the equations of motion
   real(8) Vy_e_drift      !! Dim-less electron ExB-drift y-velocity

   integer BC_flag                !! Flag, defines the boundary conditions for the potential
   integer Density_flag           !> Flag, initial distribution control:
                                  !! if = 0, uniform
                                  !! |  = 1, parabolic
                                  !! |  = 2, exponential decay
   integer N_of_particles_cell    !! Number of macroparticles of one species per cell
   integer N_of_cells_debye       !! Number of minimal cells per electron Debye length
   integer N_max_vel              !! Factor defines the maximal expected electron velocity (in thermal velocities, used for calculation of timestep)
   integer N_max_vel_distrib      !! Factor defines the maximal velocity for initial distribution (in thermal velocities)
   integer picosec_flag           !! if =/=1, then N_max_vel gives the value of time step in picoseconds (times a power of 10)
   integer micron_flag            !! if =/=0, then N_of_cells_debye gives the value of time step in microns (times a power of 10)

   integer N_box_vel              !! Number of velocity boxes per unit of initial thermal velocity
   !! Note, this value affects the scaling !
   real(8) r_max_vel, r_cells_debye !if cell size and time step are specified explicitly

   real(8) W_plasma_s1   !! The electron plasma circular frequency [s^-1]
                         !! W_plasma_s1 = SQRT(N_plasma_m3 * e_Cl**2 / (eps_0_Fm * m_e_kg))
   ! we assume that N_e = N_i and we will use only
   ! two species now - the electrons and the ions
!  real(8) W_cycl_x_s1   ! The electron cyclotron frequency for Bx [s^-1]
!  real(8) W_cycl_y_s1   ! The electron cyclotron frequency for By [s^-1]
   real(8) v_Te_ms       !! The electron thermal velocity [m/s]
                         !! v_Te_ms = SQRT(2.0_8 * T_e_eV * e_Cl / m_e_kg)
   real(8) L_debye_m     !! The Debye length [m]
                         !! L_debye_m = v_Te_ms / W_plasma_s1

   real(8) E_scl_Vm      !! Scaling factor for E [V/m], E_scl_Vm = m_e_kg * v_Te_ms * W_plasma_s1 / e_Cl
   real(8) U_scl_V       !! Scaling factor for the potential F [V]
                         !! U_scl_V = E_scl_Vm * delta_x_m

   real(8) B_scl_Tl      ! Scaling factor for B [Tesla]
   real(8) B_scl_Gs      ! Scaling factor for B [Gauss] - note, different units
   real(8) R_scl_nClm3   ! Scaling factor for charge density [nCl/m^3]
   real(8) N_scl_m3      !! Scaling factor for particles density [m^-3]
                         !! N_scl_m3 = N_plasma_m3 / DBLE(N_of_particles_cell)
   real(8) V_scl_ms      !! Scaling factor for particle velocity [m/s]
                         !! V_scl_ms = v_Te_ms / DBLE(N_box_vel)

   real(8) factor_SR     !! Dimensionless factor which appears when there is an external circuit
                         !! factor_SR = (delta_x_m * delta_t_s) / (1.0d-4 * S_electrode_cm2 * R_ext_ohm * eps_0_Fm)
                         !!           = 0.0_8 if R_ext_ohm <= 0.0_8 and S_electrode_cm2 <= 0.0_8
   real(8) factor_sigma  !! Dimensionless factor which appears when there is an external circuit
                         !! factor_sigma = 1.0_8/(e_Cl * N_plasma_m3 * delta_x_m) * dble(N_of_particles_cell)
   real(8) factor_j      !! Dimensionless factor which appears when there is an external circuit
                         !! factor_j = delta_t_s * factor_sigma

   real(8) U_ext         !! The externally maintained potential difference between the metal end plates [dim_less]
                         !! U_ext = U_ext_V / U_scl_V

   integer Q_left        !! Surface charge on the left  plasma-wall boundary [dim-less]
   integer Q_right       !! Surface charge on the right plasma-wall boundary [dim-less]
   real(8) full_Q_left   !! Surface charge on the left  plasma-wall boundary which may account for external circuit effects [dim-less]
   real(8) full_Q_right  !! Surface charge on the right plasma-wall boundary which may account for external circuit effects [dim-less]
   real(8) Q_ext         !! Surface charge deposited at x=0 ("left" )due to external current source [dim-less]

   integer N_of_part_initial   !! Initial number of macroparticles of one species

   real(8) delta_t_s     !! Timestep [s]
   real(8) delta_x_m     !! Cell size [m]

   integer T_cntr        !! Timestep counter
   integer Start_T_cntr  !! Starting timestep index
   integer Max_T_cntr    !! Number of timesteps in simulation
                         !! Max_T_cntr = 1.0d-9 * T_sim_ns / delta_t_s
   integer N_nodes       !! Number of plasma nodes (node 0 is at x = 0, and node N_nodes - 1 is at x = L_plasma_m)
   integer N_cells       !! Number of cells in plasma (N_cells = N_nodes - 1)

   integer SaveCheck_step            !! Time interval (in steps) for creating the checkpoints (full save)
   integer Restore_from_checkpoint   !! Flag, controls the initialization
   !! if 0 then - ordinary initialization of particles, diagnostics and snapshots
   !! if 1 then - particles parameters (x, v, ax, tag) are being read from the checkpoint datafiles
   !!             * assume that the diagnostic parameters for temporal output
   !!               (WriteOut_step, WriteStart_step, WriteAvg_step)
   !!               is not changed compared to the simulation, which produced the checkpoint
   !!             * assume that the parameters of snapshot creation remains unchanged
   !!               for the groups, which precede or include the checkpoint;
   !!               the subsequent groups can be modified
   !!             * the upper time limit can be modified
   !! if 2 then -  particles parameters (x, v, ax, tag) are being read from the checkpoint datafiles
   !!             * initialization of temporal diagnostics and snapshots occurs in regular manner,
   !!               i.e. on the basis of data from files "ssc_initial.dat" and "ssc_snapshot.dat" only
   character(16) check_g_filename    !! name of file, where the general data will be saved at the checkpoint
   character(16) check_e_filename    !! name of file, where the parameters of electrons will be saved at the checkpoint
   character(16) check_i_filename    !! name of file, where the parameters of ions will be saved at the checkpoint

   integer I_random_seed   !! The integer number which will be used by the random number generator through the whole program

   integer, dimension(:), allocatable :: seed !! for random_number

! SPECIES ARRAYS, HAVE FIXED DIMENSION (1:N_spec) *****
! particle parameters ------------------------

   real(8) Ms(1:N_spec)  !! normalized masses of species (by the electron mass)
   integer Qs(1:N_spec)  !! normalized charges of species (by the electron charge ABS. value)
                         !! A signed quanity
   real(8) QMs(1:N_spec) !! normalized charge-to-mass ratio of species (by that of electron ABS. value)
   real(8) VT(1:N_spec)  !! normalized thermal velocities of species (by that of the electrons)
                         !! VT(1) = N_box_vel ( = N velocity boxes per V_th_e )
                         !! VT(2) = sqrt( T_i_eV / (Ms(2) * T_e_eV) ) * N_box_vel
                         !! --"-- = ( V_th_i / V_th_e ) * N_box_vel

! precalculated coefficients -----------------

   real(8) KVx         ! Coefficients for equations of motion (X-PRE-FIN-moving)  !
   real(8) K_Q         ! Coefficients used in Poisson's equation

!  real(8) K_Xi(1:N_spec)   ! Coefficients used in interpolation of effective susceptibility
   real(8), allocatable :: K_Xi(:,:)   ! Coefficients used in interpolation of effective susceptibility

   real(8) alfa_xy(1:N_spec)   ! coefficients for calculation of alfa_x and alfa_y

!  real(8)  K11(1:N_spec) ! Coefficients (before V^{n-1/2}) for equations of motion (X-PRE-acceleration)
!  real(8)  K12(1:N_spec) ! Coefficients (before V^{n-1/2}) for equations of motion (X,Y-PRE-acceleration)
!  real(8)  K13(1:N_spec) ! Coefficients (before V^{n-1/2}) for equations of motion (X,Z-PRE-acceleration)
!  real(8)  K22(1:N_spec) ! Coefficients (before V^{n-1/2}) for equations of motion (Y-PRE-acceleration)
!  real(8)  K23(1:N_spec) ! Coefficients (before V^{n-1/2}) for equations of motion (Y,Z-PRE-acceleration)
!  real(8)  K33(1:N_spec) ! Coefficients (before V^{n-1/2}) for equations of motion (Z-PRE-acceleration)

!  real(8)  A11(1:N_spec) ! Coefficients (before Ax^{n-1}) for equations of motion (X-PRE-acceleration)
!  real(8)  A21(1:N_spec) ! Coefficients (before Ax^{n-1}) for equations of motion (Y-PRE-acceleration)
!  real(8)  A31(1:N_spec) ! Coefficients (before Ax^{n-1}) for equations of motion (Z-PRE-acceleration)

!  real(8)  A13(1:N_spec) ! Coefficients with Ez_external for equations of motion (X-PRE-acceleration)
!  real(8)  A23(1:N_spec) ! Coefficients with Ez_external for equations of motion (Y-PRE-acceleration)
!  real(8)  A33(1:N_spec) ! Coefficients with Ez_external for equations of motion (Z-PRE-acceleration)

   real(8) A_123_1             ! factor for coefficients before Ax^{n-1} for equations of motion (PRE-acceleration)
   real(8) A_123_3(1:N_spec)   ! factors for coefficients with Ez_external for equations of motion (X-PRE-acceleration)

!  real(8) KvEx(1:N_spec) ! Coefficients before Ex^{n+1} for final VX correction
!  real(8) KvEy(1:N_spec) ! Coefficients before Ex^{n+1} for final VY correction
!  real(8) KvEz(1:N_spec) ! Coefficients before Ex^{n+1} for final VZ correction

   real(8) KvE_xyz(1:N_spec)  ! factors for coefficients before Ex^{n+1} for final velocity correction

!  real(8) KxEx(1:N_spec) ! Coefficients before Ex^{n+1} for final X correction

   real(8) factor_KxEx  ! factor for coefficients before Ex^{n+1} for final X correction

! service data -------------------------------

   integer N_part_client_spec(1:256,1:N_spec)     ! Number of particles in client processes, (1-st index - process, 2-nd index - species)
   integer N_part(1:N_spec)                  !! Number of particles - total number for the server,
                                             !! number of particles in the client processes for the clients
   integer Length(1:N_spec)                  ! Lengthes of allocated arrays
   integer N_inject(1:N_spec)                ! Number (counter) of injected particles (due to ionization or SEE)

! allocatable ARRAYS *****************************
! mesh values --------------------------------
   real(8), allocatable :: EX(:)        !! The X-electric field on the mesh (0:N_cells)

   real(8), allocatable :: EX_old(:)    !! The X-electric field on the mesh (0:N_cells) at the previous time step
                                        !! This is updated the step before a snapshot is made, in order to calculate
                                        !! displacement current at the snapshot timestep

   real(8), allocatable :: GradEX(:)    !! The X-gradient of X-electric field, used for
                                        !! interpolation of EX on the mesh to particle positions
                                        !! Allocated as (0:(N_cells-1))

   real(8), allocatable :: F(:)         !! The electric potential on the mesh (0:N_cells)

   real(8), allocatable :: Xi(:)        !! The effective dielectric susceptibility, calculated from Q_strm_spec in
                                        !! CALCULATE_STR_CHARGE_DENSITY and used in CALCULATE_STR_LONG_ELECTR_FIELD

   real(8), allocatable :: Q_stream(:)  !! The instantaneous streaming charge, calculated as
                                        !! Q_stream(:) = Qs(s) * Q_strm_spec(:,s) in CALCULATE_STR_CHARGE_DENSITY
                                        !! Allocated as (0:N_cells)

   real(8), allocatable :: Q_strm_spec(:,:)  !! Instantaneous streaming charge of each species (necessary because susceptibility
                                             !! depends on mass) at node points. Accumulated in particle push procedures and
                                             !! during particle creation. Used in calculating Q_stream and Xi in 
                                             !! CALCULATE_STR_CHARGE_DENSITY, and for calculating instantaneous ion (Ni) 
                                             !! and electron (Ne) densities in diagnostics.
                                             !! Allocated as (0:N_cells,1:N_spec)
! allocatable array types --------------------

   type dalloc_array
      real(8), pointer :: part(:)
   end type dalloc_array

   type ialloc_array
      integer, pointer :: part(:)
   end type ialloc_array

!  type macroparticle
!     real(8) X
!     real(8) VX
!     real(8) VY
!     real(8) VZ
!     real(8) AX
!     integer Tag
!     real(8) old_VX
!  end type macroparticle
!
!  type macroparticle_array
!     type(macroparticle), pointer :: part(:)
!  end type macroparticle_array
!
!  type(macroparticle_array), allocatable :: species(:)
!
!### should work like this: species(s)%part(k)%X

! particles dynamics -------------------------

   type(dalloc_array), allocatable :: X_of_spec(:)       ! Array of pointers on arrays of X-coordinate of particles of blocks
   type(dalloc_array), allocatable :: VX_of_spec(:)      ! Array of pointers on arrays of X-velocity of particles of blocks
   type(dalloc_array), allocatable :: VY_of_spec(:)      ! Array of pointers on arrays of Y-velocity of particles of blocks
   type(dalloc_array), allocatable :: VZ_of_spec(:)      ! Array of pointers on arrays of Z-velocity of particles of blocks
   type(dalloc_array), allocatable :: AX_of_spec(:)      ! Array of pointers on arrays of smoothed X-acceleration of particles of blocks
   type(ialloc_array), allocatable :: Tag_of_spec(:)     !! Array of pointers on arrays of tag of particles of blocks,
                                                         !!  used to track particles which switch direction or undergo collisions
                                                         !! If = 0 : particle is "confined" within th eplasma
   type(dalloc_array), allocatable :: prev_VX_of_spec(:) ! Array of pointers on arrays of old X-velocity of particles of blocks,
                                                         !  used for specifying the tag and for calculation of the energy in diagnostics

! LINKED LIST FOR INJECTED PARTICLES *****************

   type injected_particle
      real(8) X
      real(8) VX
      real(8) VY
      real(8) VZ
      real(8) AX
      integer Tag
      type (injected_particle), pointer :: Next
   end type injected_particle

   type(injected_particle), pointer :: Inj_electron
   type(injected_particle), pointer :: Current_electron
   type(injected_particle), pointer :: Inj_ion
   type(injected_particle), pointer :: Current_ion

   integer, parameter :: eTag_Coll_Neutral       =  2      ! tag to mark electron collided with neutral atom [elastic, excitation and ionization (both)]
   integer, parameter :: eTag_Emit_Left          =  1      ! tag to mark electrons reflected / backscattered / emitted from the left wall
   integer, parameter :: eTag_Emit_Right         = -1      ! tag to mark electrons reflected / backscattered / emitted from the right wall

!  logical Tags_must_be_nullified        ! flag, is true if it is necessary to clear all tags in CLEAR_TAGS, is modified in
   ! that subroutine and in snapshot procedures

!real(8) initial_electron_energy_eV           ! for testing purposes only, initial energy of electrons, [eV]
!real(8) initial_angle_of_incidence_deg       ! for testing purposes only, initial angle of incidence, [deg]

end module CurrentProblemValues

!========================================================================================================================
! >>
!
! To add a NEW kind of collision for some species you must modify all lines containing !@#$
!
module MCCollisions

   USE CurrentProblemValues, ONLY : N_spec

! allocatable array types
   type ralloc_kind_energy
      real, pointer :: kind_energy(:,:)
   end type ralloc_kind_energy

   type ialloc_activated
      integer, pointer :: activated(:)
   end type ialloc_activated

   integer maxNcolkind_spec(1:N_spec)         !! maximal number of kinds of collisions for electrons, ions
                                              !! currently we have collisions for two plasma species only
                                              !! value is hardcoded into the code

!------- to read from file ---------
   integer Neutral_flag                   !! Flags what type of gas one is using (0 = He, 1 = Ar),
                                          !! note that Xenon is the default gas and operates
                                          !! with Ar parameters in parCollisionProc.f90
   integer Collision_flag                 !! Flags what type of collision model one is using
                                          !! 0 = Phelp's database model for Argon
                                          !! 1 = Turner Benchmark model
                                          !! anything else = original EDIPIC model

   real(8) M_neutral_amu                  !! The neutral component mass         [a.m.u.]
   real(8) N_neutral_m3                   !! The neutral component density      [m^-3]
   real(8) T_neutral_eV                   !! The neutral component temperature  [eV]

   real(8) Freq_turb_e_1_s1               !! Frequency of collisions, accounting the electron turbulence, model 1, [s^-1]
   real(8) Freq_turb_i_1_s1               !! Frequency of collisions, accounting the ion turbulence, model 1,      [s^-1]

   real(8) alpha_Vscl                     !! Scaling factor, when you take some velocity from a normalized Maxwell distribution
                                          !! and multiply the velocity by the above mentioned factor,
                                          !! you obtain the dimensionless value in electron thermal velocities

   real(8) energy_factor_eV_s(1:N_spec)   !! Scaling factor for obtaining energy of particle in electron-Volts
                                          !! from the squared particle dim-less velocity :  W_kin_eV(s) = energy_factor_eV_s(s) * v**2,
                                          !! here s - species index,
                                          !! v - dim-less absolute value of velocity (in units of v_Te / N_box_vel)

   integer Colflag_kind_spec(1:5,1:N_spec)      !! For different species [electron(*,1), ion(*,2)],                !@#$
                                                !! arrays of flags, indicating that some particular kind
                                                !! of collisions with neutrals is turned on
                                                !! 1%1 <=> e-n elastic,         model 1
                                                !! 2%1 <=> e-n excitation,      model 1
                                                !! 3%1 <=> e-n ionization,      model 1
                                                !! 4%1 <=> turbulence,          model 1
                                                !! 5%1 <=> e-n excitation,      model 2
                                                !! 1%2 <=> i-n elastic,         model 1
                                                !! 2%2 <=> i-n charge exchange, model 1
                                                !! 3%2 <=> turbulence,          model 1
                                                !! 4%2 <=> yet empty
                                                !! 5%2 <=> yet empty

   integer N_en_elast                           !! number of experimental values to be read
   real, allocatable :: Energy_en_elast_eV(:)   !! energies of known cross-secion values  for e-n elastic    collisions, [eV]
   real, allocatable :: CrSect_en_elast_m2(:)   !! experimentally measured cross-sections for e-n elastic    collisions, [m^2]

   integer N_en_excit                           !! number of experimental values to be read for excitation, model 1
   real, allocatable :: Energy_en_excit_eV(:)   !! energies of known cross-secion values  for e-n excitation collisions, [eV]
   real, allocatable :: CrSect_en_excit_m2(:)   !! experimentally measured cross-sections for e-n excitation collisions, [m^2]
   real(8) Thresh_en_excit_eV                   !! energy threshold, [eV]

   integer N_en_excit_2                         !! number of experimental values to be read for excitation, model 2
   real, allocatable :: Energy_en_excit_2_eV(:) !! energies of known cross-secion values  for e-n excitation collisions, [eV]
   real, allocatable :: CrSect_en_excit_2_m2(:) !! experimentally measured cross-sections for e-n excitation collisions, [m^2]
   real(8) Thresh_en_excit_2_eV                 !! energy threshold, [eV]

   integer N_en_ioniz                           !! number of experimental values to be read
   real, allocatable :: Energy_en_ioniz_eV(:)   !! energies of known cross-secion values  for e-n ionization collisions, [eV]
   real, allocatable :: CrSect_en_ioniz_m2(:)   !! experimentally measured cross-sections for e-n ionization collisions, [m^2]
   real(8) Thresh_en_ioniz_eV                   !! energy threshold, [eV]

   logical in_elast_flag                        !! switch, determines if i-n elastic collisions are from the default collision
                                                !! frequency=10 MHz (.true.) or from experimental cross-section data (.false.)
   integer N_in_elast                           !! number of experimental values to be read
   real, allocatable :: Energy_in_elast_eV(:)   !! energies of known cross-secion values  for i-n elastic collisions, [eV]
   real, allocatable :: CrSect_in_elast_m2(:)   !! experimentally measured cross-sections for i-n elastic collisions, [m^2]

   logical in_chrgx_flag                        !! switch, determines if i-n charge exchange collisions are from the default collision
                                                !! frequency model (.true.) or from experimental cross-section data (.false.)
   integer N_in_chrgx                           !! number of experimental values to be read
   real, allocatable :: Energy_in_chrgx_eV(:)   !! energies of known cross-secion values  for i-n charge exchange collisions, [eV]
   real, allocatable :: CrSect_in_chrgx_m2(:)   !! experimentally measured cross-sections for i-n charge exchange collisions, [m^2]


   real(8) maxEnergy_eV_spec(1:N_spec)          !! maximal value of energy of colliding plasma species, [eV]

   integer Nenergyval_spec(1:N_spec)            !! number of values of energy
                                                !! where the collisions probabilities are calculated
!------- to calculate --------------

   integer, allocatable :: Ncolkind_spec(:)           ! Number of kinds of collisions for different plasma species
   ! Ncolkind_spec(s) < maxNcolkind_spec(s)

   real(8), allocatable :: deltaEnergy_eV_spec(:)     ! energy step for different colliding plasma species, dim-less

   real, allocatable :: maxColfrac_spec(:)    ! The MAXIMAL RELATIVE number of particles of one species, !@
   ! colliding at the current time level (including NULL-collisions)

   type(ialloc_activated), allocatable :: Colkind_of_spec(:)      !! For each species, points to a 1-d array of collison indexes
                                                                  !! turned on in file "msc_partcols.dat"

   type(ralloc_kind_energy), allocatable :: Colprob_of_spec(:)    !! For different species, the pointers to the 2-d arrays of
                                                                  !! probabilities of different kinds of activated collisions
                                                                  !! versus the energy of the particle of that species

   integer, allocatable :: Numbers_of_particles(:)
   integer, allocatable :: Numbers_of_collisions(:)

!------- counters, for diagnostics -----------
   integer e_n_1_count               ! e-n, elastic, model 1
   integer e_n_2_count               ! e-n, excitation, model 1
   integer e_n_3_count               ! e-n, ionization, model 1
   integer e_t_4_count               ! e-turbulence, model 1               !@#$
   integer e_n_5_count               ! e-n, excitation, model 2
   integer i_n_1_count               ! i-n, elastic, model 1
   integer i_n_2_count               ! i-n, charge exchange, model 1
   integer i_t_3_count               ! i-turbulence, model 1               !@#$

! LINKED LIST, which stores the numbers of particles participating in collisions

   type binary_tree
      integer number
      type (binary_tree), pointer :: Larger
      type (binary_tree), pointer :: Smaller
   end type binary_tree

   type(binary_tree), pointer :: Collided_particle

end module MCCollisions

!========================================================================================================================
!
module LangevinCollisions

   integer Accounted_Langevin                         ! Flag, turns on/off the Langevin e-e collisions
   integer T_skip_Lang                                ! Timesteps to skip between the application of collisions
   integer added_cells                                ! number of cells added from each side to the cell-of-interest to
   ! in order to calculate EVDF in the cell-of-interest

!  integer N_of_left_node_Lang                        ! Left boundary for the region where teh EVDF will be calculated
!  integer N_of_right_node_Lang                       ! Right -"-
   integer W_max_Lang                                 ! Maximal velocity accounted for calculation of coefficients
   ! [in units of the scale thermal velocity]

   integer N_box_w_Lang                               ! number of velocity boxes for calculation of EVDF

   real(8) Factor_F_drag                              ! precalculated factors for the Langevin coefficients
   real(8) Factor_D_1_sqrt                            !
   real(8) Factor_D_3_sqrt                            !

   real(8), allocatable :: Fd_w(:,:)                  ! tabulated distribution function over velocity module (in the electron flow frame)
   real(8), allocatable :: F_drag(:,:)                ! tabulated drag force                                 (-"-)
   real(8), allocatable :: D_1_sqrt(:,:)              ! tabulated square root of diffusion coefficient D1 in velocity space (-"-)
   real(8), allocatable :: D_3_sqrt(:,:)              ! tabulated square root of diffusion coefficient D1 in velocity space (-"-)

   real(8), allocatable ::  V_ion_threshold(:)        ! threshold value of velocity, below which Langevin coefficients for e-i
   real(8), allocatable ::  V_ion_threshold_2(:)      ! collisions do not change
   real(8), allocatable ::  V_ion_threshold_sqrt(:)   ! is equal to the doubled initial ion thermal velocity

   real(8), allocatable :: Vx_avg(:)                  !
   real(8), allocatable :: Vy_avg(:)                  ! The averaged {flow} electron velocity
   real(8), allocatable :: Vz_avg(:)                  !
   real(8), allocatable :: energy_before(:)           ! total kinetic energy before collisions
   real(8), allocatable :: energy_after(:)            ! total kinetic energy after collisions, before the corrections are applied
   real(8), allocatable :: Log_coul(:)                ! Coulomb logarithm

   real(8), allocatable :: Density_e_sqrt(:)          ! array of square root of electron density in the middles of cells,
   ! used only by clients, the server node uses (sends to clients)
   ! the square root of the modified streaming density
end module LangevinCollisions

!========================================================================================================================
! >>
module SEEmission

   logical flag_LWsee                  !! flag, turns on/off SEE emission (both by electrons and ions) at the left wall

! the left plasma boundary can be the usual dielectric wall or the plasma source

   integer PlasmaSourceFlag            ! controls the properties of the left plasma boundary
   ! turns on / off the refluxing at the left wall (reflection with thermalization)

   integer AddInjectionFlag            ! turns on / off replenishing of particles, which were lost at the right wall, by the
   ! particles of the plasma source at the left plasma boundary

   integer const_add_N_to_inject_by_proc       ! number of macroparticles to be injected by each process aat each timestep
   ! in order to ensure additional injection with constant flux (individual for each client process)

   integer variable_add_N_to_inject    ! in order to achieve even load on all nodes,
   ! total_N_to_inject = const_add_N_to_inject_by_proc * (N_of_processes-1) + variable_add_N_to_inject
   ! the latter number is distributed between different processes all the time:
   ! one process injects one pair

! we have three possible processes when an electron hits the wall, producing a secondary electron:
! 1 - elastic backscattering
! 2 - inelastic backscattering
! 3 - true secondary emission

   integer Emitted_model(1:3)          ! for each possible process determines the way of processing

   integer Elast_refl_type             ! type of reflection for elastically reflected electron: specular / random

! ELASTIC, MODEL 1
   real(8) setD_see_elastic            ! the coefficient (ratio of emitted to incident electrons)
   real(8) minE_see_elastic_eV         ! LOWER energy boundary, [eV]
   real(8) minE_see_elastic            ! LOWER energy boundary, [dimensionless]
   real(8) maxE_see_elastic_eV         ! UPPER energy boundary, [eV]
   real(8) maxE_see_elastic            ! UPPER energy boundary, [dimensionless]

! ELASTIC, MODEL 2
   real(8) E_elast_0_eV                ! threshold energy, [eV]
   real(8) E_elast_0                   ! threshold energy, [dimensionless]
   real(8) E_elast_max_eV              ! maximum emission energy, [eV]
   real(8) E_elast_max                 ! maximum emission energy, [dimensionless]
   real(8) maxD_elast                  ! maximum emission yield
   real(8) dE_elast_eV                 ! rate of decaying (like half-width) at energies > E_elast_max_eV, [eV]
   real(8) dE_elast                    ! rate of decaying (like half-width) at energies > E_elast_max_eV, [dimensionless]
   real(8) Frac_elast_highenergy       ! fraction of total SEE yield at high energies, must be  ZERO if the model is not used

! INELASTIC, MODEL 1
   real(8) setD_see_inelastic          ! the coefficient (ratio of emitted to incident electrons)
   real(8) minE_see_inelastic_eV       ! LOWER energy boundary, [eV]
   real(8) minE_see_inelastic          ! LOWER energy boundary, [dimensionless]
   real(8) maxE_see_inelastic_eV       ! UPPER energy boundary, [eV]
   real(8) maxE_see_inelastic          ! UPPER energy boundary, [dimensionless]

! INELASTIC, MODEL 2
   real(8) Frac_inelastic              ! fraction of total SEE yield, must be ZERO if the model is not used

! TRUE SECONDARY, MODEL 1
   real(8) setD_see_true               ! the coefficient (ratio of emitted to incident electrons)
   real(8) minE_see_true_eV            ! LOWER energy boundary, [eV]
   real(8) minE_see_true               ! LOWER energy boundary, [dimensionless]
   real(8) maxE_see_true_eV            ! UPPER energy boundary, [eV]
   real(8) maxE_see_true               ! UPPER energy boundary, [dimensionless]

! TRUE SECONDARY, MODEL 2
   real(8) E_see_0_eV                  ! threshold energy, [eV]
   real(8) E_see_0                     ! threshold energy, [dimensionless]
   real(8) E_see_max_eV                ! maximal emission energy, [eV]
   real(8) E_see_max                   ! maximal emission energy, [dimensionless]
   real(8) maxD_see_classic            ! maximal emission coefficient (normal to the surface)
   real(8) k_smooth                    ! smoothness factor (from 0 for very rough to 2 for polished surfaces)

   real(8) T_see_true_eV               ! Temperature of injected true secondary electrons, [eV]

   integer electron_reflux_count       ! counter of refluxed electrons
   integer ion_reflux_count            ! counter of refluxed ions
   real(8) electron_reflux_energy      ! stores the energy of refluxed electrons
   real(8) ion_reflux_energy           ! stores the energy of refluxed ions

   integer electron_reemit_count       ! counter of refluxed electrons
   integer ion_reemit_count            ! counter of refluxed ions
   real(8) electron_reemit_energy      ! stores the energy of refluxed electrons
   real(8) ion_reemit_energy           ! stores the energy of refluxed ions

   integer N_of_lost_ions                ! counter of ions lost at the right wall (this number will be reemitted)

   integer Ion_interac_model        ! Flag (0/1), defines interaction of ions with the walls

   integer ion_left_reflect_count        ! counter of ions reflected from the wall
   real(8) ion_left_reflect_energy       ! store the energy of ions reflected from the wall

   integer ion_right_reflect_count        ! counter of ions reflected from the wall
   real(8) ion_right_reflect_energy       ! store the energy of ions reflected from the wall

   integer see_left_elastic_count        ! counter of particles, injected as elastically reflected                  tag = 1
   integer see_left_inelastic_count      ! counter of particles, injected as inelastically backscattered            tag = 1
   integer see_left_true_count           ! counter of particles, injected as true secondary                         tag = 1
   real(8) see_left_elastic_energy       ! stores the energy of particles, injected as elastically reflected
   real(8) see_left_inelastic_energy     ! stores the energy of particles, injected as inelastically backscattered
   real(8) see_left_true_energy          ! stores the energy of particles, injected as true secondary

   integer see_right_elastic_count        ! counter of particles, injected as elastically reflected                 tag = -1
   integer see_right_inelastic_count      ! counter of particles, injected as inelastically backscattered           tag = -1
   integer see_right_true_count           ! counter of particles, injected as true secondary                        tag = -1
   real(8) see_right_elastic_energy       ! stores the energy of particles, injected as elastically reflected
   real(8) see_right_inelastic_energy     ! stores the energy of particles, injected as inelastically backscattered
   real(8) see_right_true_energy          ! stores the energy of particles, injected as true secondary

   integer prie_left_from_right_count     ! counter of electrons, which come to the left wall as secondaries emitted from the right wall             tag = -1
   integer prie_left_after_coll_count     ! counter of electrons, which come to the left wall after collisions with neutrals (all kinds)             tag = 2
   real(8) prie_left_from_right_energy    ! stores the energy of electrons, which come to the left wall as secondaries emitted from the right wall
   real(8) prie_left_after_coll_energy    ! stores the energy of electrons, which come to the left wall after collisions with neutrals (all kinds)

   integer prie_right_from_left_count     ! counter of electrons, which come to the right wall as secondaries emitted from the left wall             tag = 1
   integer prie_right_after_coll_count    ! counter of electrons, which come to the right wall after collisions with neutrals (all kinds)            tag = 2
   real(8) prie_right_from_left_energy    ! stores the energy of electrons, which come to the right wall as secondaries emitted from the left wall
   real(8) prie_right_after_coll_energy   ! stores the energy of electrons, which come to the right wall after collisions with neutrals (all kinds)

   integer sece_left_from_right_count     ! counter of electrons, which are produced at the left wall by secondaries from the right wall
   integer sece_left_after_coll_count     ! counter of electrons, which are produced at the left wall by electrons after collisions with neutrals (all kinds)
   integer sece_right_from_left_count     ! counter of electrons, which are produced at the right wall by secondaries from the right wall
   integer sece_right_after_coll_count    ! counter of electrons, which are produced at the right wall by electrons after collisions with neutrals (all kinds)

end module SEEmission

!========================================================================================================================
!
module IonInducedSEEmission

! the simplest case
   real(8) setD_ionsee_true  ! constant emission coefficient, must be below 1
   real(8) minE_ionsee_eV    ! threshold energy of ions [eV], no emission if ion energy is below the threshold
   real(8) T_ionsee_eV       ! temperature of emitted electrons [eV]

! related diagnostics
   integer ionsee_left_count         ! counter of electron particles, injected at the left wall  tag = 1
   integer ionsee_right_count        ! counter of electron particles, injected at the right wall  tag = -1

   real(8) ionsee_left_energy        ! stores the energy of particles, injected at the left wall
   real(8) ionsee_right_energy       ! stores the energy of particles, injected at the right wall

end module IonInducedSEEmission

!========================================================================================================================
!
module ElectronInjection

   integer eBeamInjectFlag_left                 !! flag, turns on / off constant electron emission at the LEFT wall
   integer eBeamInjectFlag_right                !! flag, turns on / off constant electron emission at the RIGHT wall
   integer iBeamInjectFlag_left                 !! flag, turns on / off constant ion emission at the LEFT wall
   integer iBeamInjectFlag_right                !! flag, turns on / off constant ion emission at the RIGHT wall

   real(8) Delay_of_e_injection_ns_left         !!  left wall, delay between start of simulation and start of emission, [ns] ([timesteps] if negative)
   real(8) Delay_of_e_injection_ns_right        !! right wall, delay between start of simulation and start of emission, [ns] ([timesteps] if negative)
   real(8) Delay_of_i_injection_ns_left         !!  left wall, delay between start of simulation and start of emission, [ns] ([timesteps] if negative)
   real(8) Delay_of_i_injection_ns_right        !! right wall, delay between start of simulation and start of emission, [ns] ([timesteps] if negative)

   real(8) eBeam_energy_eV_left                 !!  left wall, initial energy of cold beam / temperature of warm source, [eV]
   real(8) eBeam_energy_eV_right                !! right wall, initial energy of cold beam / temperature of warm source, [eV]
   real(8) iBeam_energy_eV_left                 !!  left wall, initial energy of cold beam / temperature of warm source, [eV]
   real(8) iBeam_energy_eV_right                !! right wall, initial energy of cold beam / temperature of warm source, [eV]

   real(8) eBeam_J_Am2_left                     !!  left wall, current density of injected electron stream at the moment of injection only
   real(8) eBeam_J_Am2_right                    !! right wall, current density of injected electron stream at the moment of injection only
   real(8) iBeam_J_Am2_left                     !!  left wall, current density of injected ion stream at the moment of injection only
   real(8) iBeam_J_Am2_right                    !! right wall, current density of injected ion stream at the moment of injection only

   real(8) VX_e_beam_left                       !!  left wall, dimensionless initial velocity of electron beam
   real(8) VX_e_beam_right                      !! right wall, dimensionless initial velocity of electron beam
   real(8) VX_i_beam_left                       !!  left wall, dimensionless initial velocity of ion beam
   real(8) VX_i_beam_right                      !! right wall, dimensionless initial velocity of ion beam

   integer N_e_to_inject_total_left               !!  left wall, total number of electron macroparticles to be emitted at each timestep
   integer N_e_to_inject_total_right              !! right wall, total number of electron macroparticles to be emitted at each timestep
   integer N_i_to_inject_total_left               !!  left wall, total number of ion macroparticles to be emitted at each timestep
   integer N_i_to_inject_total_right              !! right wall, total number of ion macroparticles to be emitted at each timestep

   integer const_N_e_to_inject_by_proc_left       !!  left wall, constant number of electron particles to be emitted by each process each timestep
   integer const_N_e_to_inject_by_proc_right      !! right wall, constant number of electron particles to be emitted by each process each timestep
   integer const_N_i_to_inject_by_proc_left       !!  left wall, constant number of ion particles to be emitted by each process each timestep
   integer const_N_i_to_inject_by_proc_right      !! right wall, constant number of ion particles to be emitted by each process each timestep

   integer variable_N_e_to_inject_left            !!  left wall, variable number of electron particles to be emitted by each timestep (distributed between all processes)
   integer variable_N_e_to_inject_right           !! right wall, variable number of electron particles to be emitted by each timestep (distributed between all processes)
   integer variable_N_i_to_inject_left            !!  left wall, variable number of ion particles to be emitted by each timestep (distributed between all processes)
   integer variable_N_i_to_inject_right           !! right wall, variable number of ion particles to be emitted by each timestep (distributed between all processes)

   integer inject_e_every_this_many_timesteps_left  !!  left wall, an electron particle must be emitted every this many timesteps (1=continuous, 2=every other, etc)
   integer inject_e_every_this_many_timesteps_right !! right wall, an electron particle must be emitted every this many timesteps (1=continuous, 2=every other, etc)
   integer inject_i_every_this_many_timesteps_left  !!  left wall, an ion particle must be emitted every this many timesteps (1=continuous, 2=every other, etc)
   integer inject_i_every_this_many_timesteps_right !! right wall, an ion particle must be emitted every this many timesteps (1=continuous, 2=every other, etc)

   integer Delay_of_e_injection_T_cntr_left       !!  left wall, calculated delay between the start of simulation and the start of electron emission [timesteps]
   integer Delay_of_e_injection_T_cntr_right      !! right wall, calculated delay between the start of simulation and the start of electron emission [timesteps]
   integer Delay_of_i_injection_T_cntr_left       !!  left wall, calculated delay between the start of simulation and the start of ion emission [timesteps]
   integer Delay_of_i_injection_T_cntr_right      !! right wall, calculated delay between the start of simulation and the start of ion emission [timesteps]

   real(8) Gap_betw_e_part_left                   !!  left wall, distance between neighbour electron beam particles
   real(8) Gap_betw_e_part_right                  !! right wall, distance between neighbour electron beam particles

   integer UseSmartTagsFlag                     !! flag, controls the modification of the particle's tag related with the change of particle's velocity
                                                !! If = 1 : Particles which change direction of propagation (velocity sign) inside the plasma obtain a 
                                                !! zero tag, like the confined plasma background
                                                !! If = 0 : A change of the particle's direction of propagation will not modify the particle's tag

end module ElectronInjection

!========================================================================================================================
!
module BeamInPlasma

   integer PeriodicBoundaryFlag
   integer eBeamInPlasmaFlag          ! flag, which turns on / off the electron beam in plasma
   real(8) Start_ebeam_ns             ! delay between the start of simulation and the moment when the beam appears, [ns] ([timesteps] if negative)
   real(8) Energy_ebeam_eV            ! initial average energy of beam (energy of electron, which is at rest in the frame related with the beam), [eV]
   real(8) Te_ebeam_eV                ! temperature of the warm beam (for distribution in the frame related with the beam)
   real(8) Alfa_ebeam                 ! relative (to the scale plasma density) density of the beam
   integer N_ebeam_total              ! total number of beam macroparticles
   integer N_ebeam_proc               ! number of beam macroparticles in each process (can be different for different processes)
   integer Start_ebeam_T_cntr         ! calculated delay between the start of simulation and the appearance of the beam [timesteps]

   integer, allocatable :: N_ebeam_proc_array(:)  ! in this array the server process will store the number of beam particles in each client process
   real(8), allocatable ::  X_ebeam_part(:)       ! used during the initial distribution, X-coordinates of beam particles
   real(8), allocatable :: VX_ebeam_part(:)       ! used during the initial distribution, X-velocities of beam particles

end module BeamInPlasma

!========================================================================================================================
!
module Diagnostics
! Diagnostics follow the following format: (col. = collection)
!  __________________ : ____________________________________________________________________________________________________
!                     :
!      Procedure      :                   [   time of col.    ]         [   time of col.    ]         [   time of col.    ]
!                     : [     no col.     ]                   [ no col. ]                   [ no col. ]
!                     :                             CLEAR DATA^                   CLEAR DATA^                   CLEAR DATA^
!                     :                    AVERAGE/REPORT DATA^          AVERAGE/REPORT DATA^          AVERAGE/REPORT DATA^
!  __________________ : ____________________________________________________________________________________________________
!                     :
!      Diagnostic     :                     Finish_diag_Tcntr-v           Finish_diag_Tcntr-v           Finish_diag_Tcntr-v 
!       counters      :  Start_diag_Tcntr-v            Start_diag_Tcntr-v            Start_diag_Tcntr-v
!  __________________ : ____________________________________________________________________________________________________
!                     : -WriteStart_step--|                                                                                
!      Diagnostic     :                   |---WriteAvg_step---|         |---WriteAvg_step---|         |---WriteAvg_step---|
!      timescales     :                   |--------WriteOut_step--------|--------WriteOut_step--------|-----WriteOut_step--
!  __________________ : ____________________________________________________________________________________________________
!                     : 
!  Simulation time    : -------------------------------------------------------------------------------------------------->
!  __________________ : ____________________________________________________________________________________________________
!
! IMPORTANT: (Almost) All diagnostic data is cleared at timestep Finish_diag_Tcntr once data is written to a file. Most of 
!            this data is averaged over the length of collection (WriteAvg_step) and the averaging assumes that you collect 
!            over an interval of length WriteAvg_step. Since data is cleared at Finish_diag_Tcntr, which occurs after a 
!            period of length WriteOut_step, data can be erased if WriteAvg_step > WriteOut_step. This is not a problem if 
!            WriteAvg_step <= WriteOut_step. Graphically, the case of WriteAvg_step > WriteOut_step looks like the figure
!            below. If this is the case, diagnostics windows will overlap and lead to errors in reporting.
!
!  __________________ : ____________________________________________________________________________________________________
!                     :
!                     :                   [     time of col.   |  ]                 [     time of col.   |  ]
!      Procedure      :                                        [     time of col.   |  ]
!                     :                              CLEAR DATA^          CLEAR DATA^          CLEAR DATA^
!                     :                        AVERAGE/REPORT DATA^  AVERAGE/REPORT DATA^  AVERAGE/REPORT DATA^
!  __________________ : ____________________________________________________________________________________________________
!                     :
!      Diagnostic     :                         Finish_diag_Tcntr-v  Finish_diag_Tcntr-v  Finish_diag_Tcntr-v
!       counters      :  Start_diag_Tcntr-v   Start_diag_Tcntr-v   Start_diag_Tcntr-v   Start_diag_Tcntr-v
!  __________________ : ____________________________________________________________________________________________________
!                     : -WriteStart_step--|
!      Diagnostic     :                   |-----WriteAvg_step-----|                 |-----WriteAvg_step-----|
!      timescales     :                                        |-----WriteAvg_step-----|                 |--WriteAvg_step->
!                     :                   |----WriteOut_step---|----WriteOut_step---|----WriteOut_step---|--WriteOut_step->
!  __________________ : ____________________________________________________________________________________________________
!                     :
!  Simulation time    : -------------------------------------------------------------------------------------------------->
!  __________________ : ____________________________________________________________________________________________________

   USE CurrentProblemValues, ONLY : N_spec

   integer WriteOut_step             !! Timesteps between diagnostic ouputs
   integer WriteStart_step           !! Timestep when diagnostic data begins being collected (ie. the first Start_diag_Tcntr)
   integer WriteAvg_step             !! Number of timesteps for averaging for diagnostics and (most) snapshots
                                     !! REQUIREMENT: WriteAvg_step <= WriteOut_step
   integer TextOut_avg_prds_to_skip  ! Periods of averaging to be skipped between text outputs, 0 = each, 1 = skip one, 2 = skip two, etc.

   integer Start_diag_Tcntr          !! Timestep where the next(current) diagnostic collection window begins(began)
   integer Finish_diag_Tcntr         !! Timestep where the next(current) diagnostic collection window finishes and
                                     !! diagnostic data is cleared

   real(8) Averaging_factor          !! Inversed number of timesteps used for averaging
                                     !! Averaging_factor = 1 / WriteAvg_step
   integer N_of_saved_records        !! Counter of records in each data file (versus time)
   integer Save_check_Tcntr          !! Value of T_cntr when we create a checkpoint,
                                     !! must always coincide with a Finish_diag_Tcntr step

   integer N_of_probes                           !! number of probes for time dependencies

   integer, allocatable :: probe_node(:)         !! Locations (node numbers) of probes for time dependency diagnostics

   real(8) J_scl_Am2     !! Scaling factor for electric current density [nA/m^2]
                         !! J_scl_Am2 = Averaging_factor * (N_plasma_m3 / N_of_particles_cell) * e_Cl * (v_Te_ms / N_box_vel)

   ! This value is placed here because it includes Averaging_factor

! DIAGNOSTIC MESH-SPECIES ARRAYS ****************************

   real(8), allocatable :: VX2_cell(:,:)         !! VX^2 for each cell (space between two neighboring nodes) collcted over the
                                                 !! averaging period
   real(8), allocatable :: VY2_cell(:,:)         ! - " - VY^2
   real(8), allocatable :: VZ2_cell(:,:)         ! - " - VZ^2

   real(8), allocatable :: VX_cell(:,:)          ! - " - VX
   real(8), allocatable :: VY_cell(:,:)          ! - " - VY
   real(8), allocatable :: VZ_cell(:,:)          ! - " - VZ

   integer, allocatable :: Npart_cell(:,:)       ! - " - number of particles

   real(8), allocatable :: QF_mesh(:)            ! collects potential energy in each node (usual distribution b/w nodes)
   real(8), allocatable :: QVX_mesh(:)           ! - " - electrical current density in X-direction
   real(8), allocatable :: QVY_mesh(:)           ! - " - electrical current density in Y-direction
   real(8), allocatable :: QVZ_mesh(:)           ! - " - electrical current density in Z-direction

   real(8), allocatable :: flux_m2s1(:,:)        !! Instantaneous flux of particles at each node

!*** diagnostic arrays for ion/electron fluxes and ionization rate, 04/28/14:
   integer N_new_cell_ON_step                !! Timestep when ion-electron pair creation counting will begin.
                                             !! = 500 RF periods if RF is on, or
                                             !! = 100 ns if RF off
   integer, allocatable :: N_new_cell(:)     !! Ion-electron pair creation counter
   real(8), allocatable :: NVX_mesh(:,:)     !! Ion and electron flux counter at each node
   real(8), allocatable :: P_heat_cell(:,:)  !! Power deposited into neutral gas due to collisions with each species - 09/09/14

! energy values, dimensional
   real(8) Energy_full_eV                        ! full energy of system, [eV]

   real(8) Energy_pot_eV(1:N_spec)                    ! potential energy, [eV]
   real(8) Energy_pot(1:N_spec)                       ! - " -                                                             IS USED FOR ACCUMULATION

   real(8) Energy_kin_eV(1:N_spec)                    ! kinetic energy, 1 = electrons, 2 = ions, [eV]

   integer N_collect_part(1:N_spec)                   ! number of collected particles, 1 = electrons, 2 = ions

   real(8) Avg_kin_energy_eV(1:N_spec)                ! average kinetic energy of one particle (1=electron,2=ion), [eV]
   real(8) Avg_kin_energy_x_eV(1:N_spec)              ! - " -, X-motion,
   real(8) Avg_kin_energy_y_eV(1:N_spec)              ! - " -, Y-motion,
   real(8) Avg_kin_energy_z_eV(1:N_spec)              ! - " -, Z-motion,
   real(8) Avg_kin_energy_x(1:N_spec)                 ! - " -, X-motion,                                                 ARE USED FOR ACCUMULATION
   real(8) Avg_kin_energy_y(1:N_spec)                 ! - " -, Y-motion,                                                 ARE USED FOR ACCUMULATION
   real(8) Avg_kin_energy_z(1:N_spec)                 ! - " -, Z-motion,                                                 ARE USED FOR ACCUMULATION

   real(8) Energy_wall_eV(1:N_spec)                   ! kinetic energy of (1=electrons,2=ions), which hit the walls, [eV]
   real(8) Energy_emit_eV(1:N_spec)                   ! kinetic energy of particles emitted/reflected from the walls (1=electrons,2=ions), [eV]
   real(8) Energy_coll_eV(1:N_spec)                   ! kinetic energy of (1=electrons,2=ions), lost in collisions of electrons with neutrals, [eV]

   integer Rate_number_leftwall(1:N_spec)             ! number of (1=electrons,2=ions), which hit the  LEFT wall during the examined time interval
   integer Rate_number_rightwall(1:N_spec)            ! number of (1=electrons,2=ions), which hit the RIGHT wall during the examined time interval
   real(8) Rate_number_leftwall_ns1(1:N_spec)         ! - " -, [1/ns]
   real(8) Rate_number_rightwall_ns1(1:N_spec)        ! - " -, [1/ns]

   integer Rate_number_leftemit(1:N_spec)             ! number of (1=electrons,2=ions), which were emitted from the  LEFT wall during the examined time interval
   integer Rate_number_rightemit(1:N_spec)            ! number of (1=electrons,2=ions), which were emitted from the RIGHT wall during the examined time interval
   real(8) Rate_number_leftemit_ns1(1:N_spec)         ! - " -, [1/ns]
   real(8) Rate_number_rightemit_ns1(1:N_spec)        ! - " -, [1/ns]

   real(8) Yield_left(1:N_spec)                       ! relative emission coefficient (1=electrons,2=ions), LEFT wall
   real(8) Yield_right(1:N_spec)                      ! relative emission coefficient (1=electrons,2=ions), RIGHT wall
   real(8) eYield_left_elast,  eYield_left_inelast,  eYield_left_true   ! relative coefficients for components of electron emission
   real(8) eYield_right_elast, eYield_right_inelast, eYield_right_true  !

   real(8) Rate_energy_leftwall(1:N_spec)                  ! kinetic energy of (1=electrons,2=ions), which hit the  LEFT wall
   real(8) Rate_energy_rightwall(1:N_spec)                 ! kinetic energy of (1=electrons,2=ions), which hit the RIGHT wall
   real(8) Avg_energy_leftwall_eV(1:N_spec)           ! average kinetic energy of (1=electrons,2=ions), which hit the  LEFT wall, [eV]
   real(8) Avg_energy_rightwall_eV(1:N_spec)          ! average kinetic energy of (1=electrons,2=ions), which hit the RIGHT wall, [eV]

   real(8) Rate_energy_leftemit(1:N_spec)                  ! kinetic energy of (1=electrons,2=ions), which were emitted from the  LEFT wall during the time interval
   real(8) Rate_energy_rightemit(1:N_spec)                 ! kinetic energy of (1=electrons,2=ions), which were emitted from the RIGHT wall during the time interval
   real(8) Avg_energy_leftemit_eV(1:N_spec)           ! average kinetic energy - " -, [eV]
   real(8) Avg_energy_rightemit_eV(1:N_spec)          ! average kinetic energy - " -, [eV]

   real(8) Energy_heat_eV(1:N_spec)                        ! energy, transferred into the plasma due to Joule heating (electric current), [eV]

   real(8) Rate_energy_full_eVns1                ! rate of change of full energy, [eV/ns]
   real(8) Rate_energy_pot_eVns1(1:N_spec)            ! rate of change of potential energy, [eV/ns]
   real(8) Rate_energy_kin_eVns1(1:N_spec)            ! rate of cange of kinetic energy, 1 = electrons, 2 = ions, [eV/ns]

   real(8) Rate_energy_wall_eVns1(1:N_spec)           ! rate of change of (1=electrons,2=ions) energy due to collisions with the walls, [eV/ns]
   real(8) Rate_energy_emit_eVns1(1:N_spec)           ! rate of change of energy due to emission of particles from the wall (1=electrons,2=ions), [eV/ns]
   real(8) Rate_energy_coll_eVns1(1:N_spec)           ! rate of change of (1=electrons,2=ions) energy due to collisions with neutral particles, [eV/ns]

   real(8) Rate_energy_heat_eVns1(1:N_spec)           ! rate of change of energy due to Joule heating, [eV/ns]

   real(8) U_app_V, U_cap_V, Sigma_electrode, Enr_consumed_Jm2 !** to calculate power for capacitive discharge

! energy rates, dim-less
   real(8) Rate_energy_wall(1:N_spec)                 ! rate of change of energy of (1=electrons,2=ions) due to collisions with the walls
   real(8) Rate_energy_emit(1:N_spec)                 ! rate of change of energy of (1=electrons,2=ions) emitted/reflected from the walls
   real(8) Rate_energy_coll(1:N_spec)                 ! rate of change of energy of (1=electrons,2=ions) due to collisions with neutral particles
   real(8) Rate_energy_heat(1:N_spec)                      ! rate of change of energy due to Joule heating

! flow values
   real(8) VY_e_avg_ms                           ! averaged (over all particles) Y-velocity for electrons, [m/s]           ARE USED FOR ACCUMULATION
   real(8) VZ_e_avg_ms                           ! averaged (over all particles) Z-velocity for electrons, [m/s]           ARE USED FOR ACCUMULATION
   real(8) JY_avg_Am2                            ! averaged (over cross section) Y-electric current density, [A/m^2]
   real(8) JZ_avg_Am2                            ! averaged (over cross section) Z-electric current density, [A/m^2]

! density
   real(8) N_avg_m3(1:N_spec)                         ! averaged (over cross section) electron (1) and ion(2) density, [1/m^3]

   real(8) prev_Q_left   ! Surface charge on the left  plasma-wall boundary at previous diagnostics time-step [dim-less]
   real(8) prev_Q_right  ! Surface charge on the right plasma-wall boundary at previous diagnostics time-step [dim-less]

! dim-less quantities accumulated only during one timestep (the last timestep of the averaging period)
   real(8) VY_recent(1:N_spec)                        ! for Y-velocity, separate for electrons (1) and ions (2)
   real(8) VZ_recent(1:N_spec)                        ! for Z-velocity, separate for electrons (1) and ions (2)

   ! energy factors [in units of 1/m^2]
   real(8) N_in_macro                           !! number of electrons/ions in one macroparticle [1/m^2]
                                                !! N_in_macro = N_plasma_m3 * delta_x_m / N_of_particles_cell
   real(8) Factor_energy_eV                     !! Factor_energy_eV = T_e_eV / N_box_vel**2
   real(8) Factor_energy_macro_eV               !! Factor_energy_macro_eV = N_in_macro * Factor_energy_eV [eV/m^2]
   real(8) Factor_energy_pot_eV                 !! Factor_energy_pot_eV = N_in_macro * (2.0_8 * T_e_eV / r_cells_debye) * 0.5_8
   real(8) Factor_rate_ns1                      !! Factor_rate_ns1 = 1.0_8 / (WriteOut_step * delta_t_s * 1.0e9) [1/ns]
   real(8) Factor_rate_macro_ns1                !! Factor_rate_macro_ns1 = N_in_macro * Factor_rate_ns1 [1/(ns*m^2)]
   real(8) Factor_rate_eVns1                    !! Factor_rate_eVns1 = Factor_energy_eV * Factor_rate_macro_ns1 [eV/(ns*m^2)]
   real(8) Factor_Joule_energy_eV               !! Factor_Joule_energy_eV = N_in_macro * (v_Te_ms / N_box_vel) * E_z_ext_Vm * (delta_t_s)

   real(8) Init_energy_full_eV

   integer text_output_counter                   ! is used for skipping text output

   real(8) f_factor(1:N_spec)                   !! Factor for calculation of collision frequencies between electrons/ions and neutrals
                                                !! f_factor(s) = 1.0_8 / (N_part(s) * WriteOut_step * delta_t_s)

   integer, parameter :: length_of_fbufer=600    !! length of the fmid and fwall buffers
   real(8) fmid_bufer(1:length_of_fbufer)        !! bufer for calculation of time-averaged potential at the midplane
   real(8) fwall_bufer(1:length_of_fbufer)       !! bufer for calculation of time-averaged potential at the left wall (x=0)

end module Diagnostics

!========================================================================================================================
!
module Snapshots
! Snapshot timing note:
!  - Snapshot data is outputted on Finish_diag_Tcntr steps, before diagnostic data (much of which is also used for 
!    snapshots) is cleared. This being the case, as long as WriteAvg_step <= WriteOut_step as described in the 
!    Diagnostics module, the snapshot data will be reported correctly.
!
   integer N_of_all_snaps                        !! Number of snapshots (counted up from zero in CALCULATE_SNAPSHOT_TIMES)
   integer current_snap                          ! index of current snapshot (which must be created)
   integer counter_of_profiles_to_save           ! number of profiles of the potential to be saved
   character(4) snapnumber_txt                   ! prefix in the snapshot files, describes the number of the snapshot

   integer Flag_pp_creation                      ! flag, turns on/off the creation of phase planes using the fraction of all particles
   integer N_to_skip                             ! the fraction mention above is controlled by this parameter (number of particles to skip) for bulk electrons
   integer N_to_skip_left                        ! -""- for electrons emitted at the left wall
   integer N_to_skip_right                       ! -""- for electrons emitted at the right wall
   integer N_to_skip_ion                         ! -""- for ions

   integer N_of_all_vdf_locs                     ! number of locations for calculations of velocity distribution functions (VDF)
   integer N_of_all_edf_locs                     ! number of locations for calculations of energy distribution functions (EDF)

   real(8) Ve_x_max               ! Maximal x-velocity (normal) for electron velocity distribution (initially given in V_therm)
   real(8) Ve_yz_max              ! Maximal y,z-velocity (parallel) for electron velocity distribution (initially given in V_therm)

   real(8) Ee_max_eV              !! Maximal energy for electron energy distribution (initially given in eV)
   real(8) Ei_max_eV              !! Maximal energy for ion energy distribution (initially given in eV)
   real(8) Ei_wall_max_eV         !! Maximal energy for ion energy distribution at the wall (initially given in eV)
   integer N_E_bins               !! Number of energy bins for electron and ion energy distributions (bin indexes run 1:N_E_bins)

   real(8) delta_Ee_eV            ! Energy step for electron energy distribution (initially given in eV)
   real(8) delta_Ei_eV            ! Energy step for ion energy distribution (initially given in eV)
   real(8) delta_Ei_wall_eV       ! Energy step for ion energy distribution at the wall (initially given in eV)

   integer N_box_Vx_e             ! (Total number of velocity boxes - 1) for positive x-velocity (normal) for electrons
   integer N_box_Vyz_e            ! (Total number of velocity boxes - 1) for positive y,z-velocity (parallel) for electrons
   integer N_box_Vx_e_low         ! = -N_box_Vx_e                   !
   integer N_box_Vx_e_top         ! =  N_box_Vx_e + 1               ! this will save three operations of addition and
   integer N_box_Vyz_e_low        ! = -N_box_Vyz_e                  ! three operations of sign changing for each electron
   integer N_box_Vyz_e_top        ! =  N_box_Vyz_e + 1              !

   integer N_box_Vx_i             ! (Number of velocity boxes - 1) for positive x-velocity (normal) for ions
   integer N_box_Vx_i_low         ! = -N_box_Vx_i                   ! this will save one addition and one sign changing for each ion
   integer N_box_Vx_i_top         ! =  N_box_Vx_i + 1

   integer, allocatable :: Tcntr_snapshot(:)     ! timesteps when the snapshot files are written

   integer, allocatable :: Vdf_location_bnd(:)   ! Spatial boundaries, defining the locations for calculation of VDFs (numbers of nodes)
   integer, allocatable :: Edf_location_bnd(:)   ! Spatial boundaries, defining the locations for calculation of EDFs (numbers of nodes)

   logical Accumulate_wall_df       !! flag, showing whether it is necessary (or not) to add particles to the wall distribution functions
   logical flag_evxdf               !! flag, controls snapshot distribution function creation
   logical flag_evydf               !! flag, controls snapshot distribution function creation
   logical flag_evzdf               !! flag, controls snapshot distribution function creation
   logical flag_ivxdf               !! flag, controls snapshot distribution function creation
   logical flag_eedf                !! flag, controls snapshot distribution function creation
   logical flag_iedf                !! flag, controls snapshot distribution function creation
   logical flag_ilwedf              !! flag, controls snapshot distribution function creation
   logical flag_irwedf              !! flag, controls snapshot distribution function creation

! DIAGNOSTIC DISTRIBUTION FUNCTION ARRAYS *******************

   real(8), allocatable :: evx_mid_of_box(:)     ! middles of   x-velocity boxes (in units of V_th_e), electrons
   real(8), allocatable :: evyz_mid_of_box(:)    ! middles of y,z-velocity boxes (in units of V_th_e), electrons

   real(8), allocatable :: eedf_mid_of_box_eV(:)      ! middles of energy boxes (in units of eV), electrons
   real(8), allocatable :: iedf_mid_of_box_eV(:)      ! middles of energy boxes (in units of eV), ions
   real(8), allocatable :: iedf_wall_mid_of_box_eV(:) ! middles of energy boxes (in units of eV), ions at the wall

   integer, allocatable :: ep_2vdf_lw(:,:)       ! 2-d (|v_norm|,|v_par|) distribution of incident (primary) electrons at the  left wall
   integer, allocatable :: ep_2vdf_rw(:,:)       ! 2-d (|v_norm|,|v_par|) distribution of incident (primary) electrons at the right wall
   integer, allocatable :: es_2vdf_lw(:,:)       ! 2-d (|v_norm|,|v_par|) distribution of emitted electrons at the  left wall
   integer, allocatable :: es_2vdf_rw(:,:)       ! 2-d (|v_norm|,|v_par|) distribution of emitted (secondary) electrons at the right wall

! common
!  integer, allocatable :: e_2vxvydf_loc(:,:,:)    ! 2-d (v_x,v_y) distribution of electrons at certain location (region) inside the plasma
!  integer, allocatable :: e_2vxvzdf_loc(:,:,:)    ! 2-d (v_x,v_z) distribution of electrons at certain location (region) inside the plasma
   integer, allocatable :: e_vxdf_loc(:,:)         !! distribution of electrons over v_x at certain location inside the plasma
   integer, allocatable :: e_vydf_loc(:,:)         !! distribution of electrons over v_y at certain location inside the plasma
   integer, allocatable :: e_vzdf_loc(:,:)         !! distribution of electrons over v_z at certain location inside the plasma
! emitted from the left wall, which do not mix with plasma
!  integer, allocatable :: ebl_2vxvydf_loc(:,:,:)  ! 2-d (v_x,v_y) distribution of electrons at certain location (region) inside the plasma
!  integer, allocatable :: ebl_2vxvzdf_loc(:,:,:)  ! 2-d (v_x,v_z) distribution of electrons at certain location (region) inside the plasma
   integer, allocatable :: ebl_vxdf_loc(:,:)       !! distribution of tagged electrons over v_x at certain location inside the plasma
   integer, allocatable :: ebl_vydf_loc(:,:)       !! distribution of tagged electrons over v_y at certain location inside the plasma
   integer, allocatable :: ebl_vzdf_loc(:,:)       !! distribution of tagged electrons over v_z at certain location inside the plasma
! emitted from the right wall, which do not mix with plasma
!  integer, allocatable :: ebr_2vxvydf_loc(:,:,:)  ! 2-d (v_x,v_y) distribution of electrons at certain location (region) inside the plasma
!  integer, allocatable :: ebr_2vxvzdf_loc(:,:,:)  ! 2-d (v_x,v_z) distribution of electrons at certain location (region) inside the plasma
   integer, allocatable :: ebr_vxdf_loc(:,:)       !! distribution of tagged electrons over v_x at certain location inside the plasma
   integer, allocatable :: ebr_vydf_loc(:,:)       !! distribution of tagged electrons over v_y at certain location inside the plasma
   integer, allocatable :: ebr_vzdf_loc(:,:)       !! distribution of tagged electrons over v_z at certain location inside the plasma

!*** diagnostic array for ion distribution at the right wall, 08/01/14:
   real(8), allocatable :: irwedf(:), ilwedf(:)  ! ion distribution colliding with the right and left walls, respectively

   real(8), allocatable :: ivx_mid_of_box(:)     !! middles of x-velocity boxes (in units of V_th_e), ions
   integer, allocatable :: i_vxdf_loc(:,:)       !! distribution of ions over v_x at certain location inside the plasma
   integer, allocatable :: i_vydf_loc(:,:)       !! distribution of ions over v_y at certain location inside the plasma
   integer, allocatable :: i_vzdf_loc(:,:)       !! distribution of ions over v_z at certain location inside the plasma

   integer, allocatable :: e_edf_loc(:,:)        !! distribution of ions over energy at certain location inside the plasma
   integer, allocatable :: i_edf_loc(:,:)        !! distribution of ions over energy at certain location inside the plasma

end module Snapshots

!----------------------------------
module TestParticles

   integer Leftmost_initial_test_part_node       ! left boundary of the initial distribution of test particles
   integer N_of_initial_test_part_spatial_layers ! number of spatial layers in the initial distribution
   real(8) d_spatial                             ! difference between initial positions of neighbor layers

   real(8) Min_initial_test_part_energy_eV       ! minimal initial energy of test particles [eV]
   integer N_of_initial_test_part_energy_layers  ! number of energy layers in the initial distribution
   real(8) d_energy_eV                           ! difference between initial energies of neighbor layers [eV]

   integer N_of_test_parts                       ! number of test particles
   integer N_of_test_part_runs                   ! number of times when the test particles will be launched
   integer current_test_run                      ! the number of the current test run

   integer N_of_test_part_save_skip              ! number of timesteps to be skipped between saving test particles
   integer test_part_skip_save_counter           ! counter of timesteps when saving the test particles was skipped

   type TestParticle
      real(8) x
      real(8) vx
      real(8) ax
      integer active              !  1 = active, 0 = non-active
   end type TestParticle

! allocatable arrays
   integer, allocatable :: Test_part_run_limits(:,:)    ! initial and final time steps for each run

   type(TestParticle), allocatable :: Test_part(:)            ! test particle parameters

end module TestParticles
