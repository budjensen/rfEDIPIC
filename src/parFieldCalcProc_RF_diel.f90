!============================================================================
! This subroutine calculates the mesh values of the streaming charge density
!
SUBROUTINE CALCULATE_STR_CHARGE_DENSITY

   use mpi
   USE ParallelOperationValues
   USE CurrentProblemValues
   IMPLICIT NONE

!  INCLUDE 'mpif.h'

!  INTEGER s                             ! species
   INTEGER ierr, ALLOC_ERR, DEALLOC_ERR
   REAL(8), ALLOCATABLE :: rbufer(:)     ! for Q_strm_spec
   REAL(8), ALLOCATABLE :: rbufer2(:)
   real(8) j_Am2_vst
   INTEGER ibufer(1:2)                   ! for Q_left and Q_right
   INTEGER ibufer2(1:2)

   INTEGER i

!  LOGICAL flag
!  INTEGER stattus(MPI_STATUS_SIZE)
!  INTEGER n, source


   IF (.NOT.ALLOCATED(rbufer)) THEN
      ALLOCATE(rbufer(1:(2*N_nodes)), STAT = ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : CALCULATE_STR_CHARGE_DENSITY : Error in ALLOCATE rbufer !!!")', Rank_of_process
         PRINT '(2x,"The program will be terminated now :(")'
         STOP
      END IF
   END IF

   IF (.NOT.ALLOCATED(rbufer2)) THEN
      ALLOCATE(rbufer2(1:(2*N_nodes)), STAT = ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : CALCULATE_STR_CHARGE_DENSITY : Error in ALLOCATE rbufer2 !!!")', Rank_of_process
         PRINT '(2x,"The program will be terminated now :(")'
         STOP
      END IF
   END IF

   rbufer2 = 0.0_8
   ibufer2 = 0

   IF (Rank_of_process.GT.0) THEN                                                             ! client >>>

! send Q_left and Q_right to the server process
      ibufer(1) = Q_left
      ibufer(2) = Q_right
      CALL MPI_REDUCE(ibufer, ibufer2, 2, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

! send Q_strm_spec to the server
      rbufer(1:N_nodes) = Q_strm_spec(0:N_cells, 1)
      IF (N_spec.EQ.1) THEN
         CALL MPI_REDUCE(rbufer, rbufer2, N_nodes, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      ELSE
         rbufer((N_nodes+1):(2*N_nodes)) = Q_strm_spec(0:N_cells, 2)
         CALL MPI_REDUCE(rbufer, rbufer2, 2*N_nodes, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      END IF
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

   ELSE                                                                                       ! server >>>

! receive (accumulate) Q_left and Q_right from all client processes
      ibufer = 0
      CALL MPI_REDUCE(ibufer2, ibufer, 2, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      Q_left  = ibufer(1)
      Q_right = ibufer(2)

      if (BC_flag.EQ.3) then !*** total time-integrated charge passed through the circuit
         Q_ext = Q_ext + factor_j * j_Am2_vst(T_cntr * delta_t_s)
      else
         Q_ext = 0.0_8
      end if

      IF (BC_flag.NE.0 .and. BC_flag.NE.2) THEN
! floating wall or specified charge moved through or U_ext behind dielectric layer, direct accumulation
         full_Q_left  = DBLE(Q_left)  + Q_ext !** continously integrated from the initial moment
         full_Q_right = DBLE(Q_right) - Q_ext
      ELSE
! constant given potential or external circuit, the accumulated value is the correction
         full_Q_left  = full_Q_left  + DBLE(Q_left)
         full_Q_right = full_Q_right + DBLE(Q_right)
      END IF

      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

! receive (accumulate) Q_strm_spec from all client processes
      rbufer = 0.0_8
      IF (N_spec.EQ.1) THEN
         CALL MPI_REDUCE(rbufer2, rbufer, N_nodes, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
         Q_strm_spec(0:N_cells, 1) = rbufer(1:N_nodes)
         Q_stream(0:N_cells) = Qs(1) * Q_strm_spec(0:N_cells, 1)
! ###
         DO i = 0, N_cells
            Xi(i) = K_Xi(i,1) * Q_strm_spec(i,1)
         END DO
      ELSE
         CALL MPI_REDUCE(rbufer2, rbufer, 2*N_nodes, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
         Q_strm_spec(0:N_cells, 1) = rbufer(1:N_nodes)
         Q_strm_spec(0:N_cells, 2) = rbufer((N_nodes+1):(2*N_nodes))
         Q_stream(0:N_cells) = Qs(1) * Q_strm_spec(0:N_cells, 1) + Qs(2) * Q_strm_spec(0:N_cells, 2)
! ###
         DO i = 0, N_cells
            Xi(i) =   K_Xi(i,1) * Q_strm_spec(i,1) + K_Xi(i,2) * Q_strm_spec(i,2)
         END DO
      END IF

   END IF

   IF (ALLOCATED(rbufer)) THEN
      DEALLOCATE(rbufer, STAT = DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : CALCULATE_STR_CHARGE_DENSITY : Error in DEALLOCATE rbufer !!!")', Rank_of_process
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

   IF (ALLOCATED(rbufer2)) THEN
      DEALLOCATE(rbufer2, STAT = DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : CALCULATE_STR_CHARGE_DENSITY : Error in DEALLOCATE rbufer2 !!!")', Rank_of_process
         PRINT *, 'The program will be terminated now :('
         STOP
      END IF
   END IF

END SUBROUTINE CALCULATE_STR_CHARGE_DENSITY

!=====================================================================================
! This subroutine calculates the mesh values of the potential F(1:N_nodes)
! at the integer nodes. Then it calculates the mesh values of the longitudinal electric field
! EZ(1:N_nodes)
! Besides, the energy of the longitudinal electric field is calculated too. IN FUTURE!!!
SUBROUTINE CALCULATE_STR_LONG_ELECTR_FIELD

   use mpi
   USE ParallelOperationValues
   USE CurrentProblemValues
   USE BeamInPlasma, ONLY : PeriodicBoundaryFlag
   USE Diagnostics, ONLY: N_in_macro
   IMPLICIT NONE

!  INCLUDE 'mpif.h'

!  INTEGER igrid

   REAL(8) a_eq(0:N_cells,1:4)
   REAL(8) b_eq(0:N_cells)

   INTEGER i, j
   REAL(8) tmp
   INTEGER ierr

! function
   REAL(8) U_ext_vst, j_Am2_vst, sigma_Cm2_vst, sawtooth_Cm2_vst, QDC_Cm2_vst

!  K_Q            = 1.0_8 / REAL(N_of_particles_cell * N_of_cells_debye), set elsewhwre
!  factor_SR = (delta_x_m * delta_t_s) / (1.0d-4 * S_electrode_cm2 * R_ext_ohm * eps_0_Fm), set elsewhere
!  factor_j = delta_t_s / (e_Cl * N_plasma_m3 * delta_x_m) * dble(N_of_particles_cell)
!  factor_sigma = 1.0_8/(e_Cl * N_plasma_m3 * delta_x_m) * dble(N_of_particles_cell), set elsewhere
!  factor_j = delta_t_s * factor_sigma
!  factor_j = delta_t_s / (e_Cl * N_in_macro)
!*** note that factor_SR / U_scl_V * (1.0d-4 * S_electrode_cm2 * R_ext_ohm) = factor_j * K_Q

   IF (Rank_of_process.EQ.0) THEN                                                          ! server >>>

      SELECT CASE (BC_flag)
       CASE (0)              ! potentials of walls (plasma boundaries) are given
         a_eq( 0, 1) =  0.0_8
         a_eq( 0, 2) =  1.0_8
         a_eq( 0, 3) =  0.0_8
         b_eq( 0)    =  U_ext_vst(T_cntr * delta_t_s)

       CASE (1)              ! potentials of walls (plasma boundaries) are floating
         a_eq( 0, 1) =  0.0_8
         a_eq( 0, 2) =  1.0_8
         a_eq( 0, 3) = -1.0_8
         b_eq( 0)    = K_Q * (full_Q_left + Q_stream(0)) / (1.0_8 + Xi(0) + Xi(1))

       CASE (2)              ! there is an external contour with a battery and a resistor
         a_eq( 0, 1) =  0.0_8
         a_eq( 0, 2) =  1.0_8 + Xi(0) + Xi(1) + factor_SR * (1.0_8 + 2.0_8 * Xi(0))
         a_eq( 0, 3) = -1.0_8 - Xi(0) - Xi(1)
         b_eq( 0)    = factor_SR * U_ext_vst(T_cntr * delta_t_s) * (1.0_8 + 2.0_8 * Xi(0)) + K_Q * (full_Q_left + Q_stream(0))

       CASE (3)              ! electrode current is specified
!!        full_Q_left = full_Q_left + factor_j * j_Am2_vst(T_cntr * delta_t_s) *** obsolete, 9/1/15
         a_eq( 0, 1) =  0.0_8
         a_eq( 0, 2) =  1.0_8
         a_eq( 0, 3) = -1.0_8
!!        b_eq( 0)    =  K_Q * ( factor_sigma * sawtooth_Cm2_vst(T_cntr * delta_t_s) + full_Q_left + Q_stream(0)) / (1.0_8 + Xi(0) + Xi(1))
         b_eq( 0)    =  K_Q * (full_Q_left + Q_stream(0)) / (1.0_8 + Xi(0) + Xi(1))

       CASE (4)              !*** dielectric layer with eps_diel at -d < x < 0, external potential specified at x = -d
         a_eq( 0, 1) =  0.0_8
         a_eq( 0, 2) =  1.0_8 + Xi(0) + Xi(1) + eps_param
         a_eq( 0, 3) = -1.0_8 - Xi(0) - Xi(1)
         b_eq( 0)    = K_Q * (full_Q_left + Q_stream(0)) + eps_param * U_ext_vst(T_cntr * delta_t_s)
!!        b_eq( 0)    = K_Q * (dble(Q_left) + Q_stream(0)) + eps_param * U_ext_vst(T_cntr * delta_t_s)
      END SELECT

      DO i = 1, N_cells - 1
         a_eq(i, 1) = 1.0_8 + Xi(i-1) + Xi(i)
         a_eq(i, 3) = 1.0_8 + Xi(i+1) + Xi(i)
         a_eq(i, 2) = - a_eq(i, 1) - a_eq(i, 3)
         b_eq(i)    = - K_Q * Q_stream(i)
      END DO

      a_eq(N_cells, 1) =  0.0_8
      a_eq(N_cells, 2) =  1.0_8
      a_eq(N_cells, 3) =  0.0_8
      b_eq(N_cells) = 0.0_8

      a_eq(0:N_cells, 4) = 0.0_8

! diagonalization . . .
      DO i = 1, N_cells !N_nodes
         IF (a_eq(i-1,2).NE.0.0_8) THEN
            tmp = a_eq(i, 1) / a_eq(i-1, 2)
            a_eq(i, 1) = 0.0_8
            a_eq(i, 2) = a_eq(i, 2) - tmp * a_eq(i-1, 3)
            b_eq(i)    = b_eq(i)    - tmp * b_eq(i-1)
         ELSE
            IF (a_eq(i-1,3).EQ.0.0_8) THEN
               PRINT *, 'Error in Poisson system!', T_cntr
               PRINT '(2x,"Equation ",i5,": a(1) =",e12.5," a(2) =",e12.5," a(3) =",e12.5," b = ",e12.5)', &
               & i, a_eq(i,1), a_eq(i,2), a_eq(i,3), b_eq(i)
               OPEN (41, FILE = 'error_a_b_eq.dat')
               DO j = -1, N_nodes
                  WRITE (41, '(2x,i5,4(2x,e12.5),6x,e12.5)') j, a_eq(j,1), a_eq(j,2), a_eq(j,3), a_eq(j,4), b_eq(j)
               END DO
               CLOSE (41, STATUS = 'KEEP')
               PRINT *, 'Program will be terminated now :('
               STOP
            END IF
            a_eq(i-1, 2) = a_eq(i, 1)
            a_eq(i, 1)   = 0.0_8

            tmp          = a_eq(i-1, 3)
            a_eq(i-1, 3) = a_eq(i, 2)
            a_eq(i, 2)   = tmp

            a_eq(i-1, 4) = a_eq(i, 3)
            a_eq(i, 3)   = 0.0_8

            tmp       = b_eq(i-1)
            b_eq(i-1) = b_eq(i)
            b_eq(i)   = tmp
         END IF
      END DO

      F(N_cells)   =  b_eq(N_cells) / a_eq(N_cells, 2)
      F(N_cells-1) = (b_eq(N_cells-1) - a_eq(N_cells-1, 3) * F(N_cells)) / a_eq(N_cells-1, 2)

      DO i = N_cells - 2, 0, -1
         F(i) = (b_eq(i) - a_eq(i, 3) * F(i+1) - a_eq(i, 4) * F(i+2)) / a_eq(i, 2)
      END DO

      EX(0)       = (-0.5_8 * K_Q * Q_stream(0)       + (F(0)         - F(1)      ) * (1.0_8 + Xi(0)         + Xi(1))       ) / (1.0_8 + 2.0_8 * Xi(0))
      EX(N_cells) = ( 0.5_8 * K_Q * Q_stream(N_cells) + (F(N_cells-1) - F(N_cells)) * (1.0_8 + Xi(N_cells-1) + Xi(N_cells)) ) / (1.0_8 + 2.0_8 * Xi(N_cells))

      DO i = 1, N_cells - 1
         EX(i) = 0.5_8 * (F(i - 1) - F(i + 1))
      END DO

      IF (PeriodicBoundaryFlag.EQ.1) THEN
         EX(0) = 0.5_8 * (F(N_cells-1) - F(1))
         EX(N_cells) = EX(0)
      END IF

   END IF

   CALL MPI_BCAST(EX, (N_cells+1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)      ! server sends and clients receive array EX

   DO i = 0, N_cells - 1
      GradEX(i) = EX(i+1) - EX(i)
   END DO

END SUBROUTINE CALCULATE_STR_LONG_ELECTR_FIELD

!-----------------------------------------------
!
REAL(8) FUNCTION U_ext_vst(t_s)

   USE CurrentProblemValues, ONLY : U_ext, U_scl_V, rf_on, U_rf_V, f_rf_Hz, t_start_s

   IMPLICIT NONE

   REAL(8) t_s

   IF (rf_on) THEN

      U_ext_vst = U_ext

      IF (t_s.LT.t_start_s) RETURN

      U_ext_vst = U_ext + (U_rf_V / U_scl_V) * SIN(6.283185307_8 * (t_s - t_start_s) * f_rf_Hz)

   ELSE
      U_ext_vst = U_ext

   END IF

END FUNCTION U_ext_vst

!-----------------------------------------------
!
!REAL(8) FUNCTION j_Am2_vst(t_s)
!
!  IMPLICIT NONE
!
!  REAL(8), PARAMETER :: j_rf_Am2 = 5.d3
!  REAL(8), PARAMETER :: f_rf_Hz = 2.712d7
!  REAL(8), PARAMETER :: t_start_s = 0.0_8
!
!  REAL(8) t_s
!
!  j_Am2_vst = 0.0_8
!
!  IF (t_s.LT.t_start_s) RETURN
!
!  j_Am2_vst =  j_rf_Am2 * cos(6.283185307_8 * (t_s - t_start_s) * f_rf_Hz)
!
!END FUNCTION j_Am2_vst

REAL(8) FUNCTION sigma_Cm2_vst(t_s)

   IMPLICIT NONE

   REAL(8), PARAMETER :: j_max_Am2 = 100.0_8
   REAL(8), PARAMETER :: f_rf_Hz = 2.712d7
   REAL(8), PARAMETER :: t_start_s = 0.0_8
   real(8) omega
   REAL(8) t_s

   omega = 6.283185307_8 * f_rf_Hz
   sigma_Cm2_vst = 0.0_8

   IF (t_s.LT.t_start_s) RETURN

   sigma_Cm2_vst =  (j_max_Am2 / omega) * SIN( omega * (t_s - t_start_s) )

END FUNCTION sigma_Cm2_vst

REAL(8) FUNCTION sawtooth_Cm2_vst(t_s)

   IMPLICIT NONE
   REAL(8), PARAMETER :: j_max_Am2 = 200.0_8
   REAL(8), PARAMETER :: f_rf_Hz = 2.712d7
   REAL(8), PARAMETER :: t_start_s = 0.0_8
   REAL(8), PARAMETER :: asymm = 0.7_8
   real(8) omega, phase, pi, twopi
   REAL(8) t_s, T_rf_s, arg, c
   integer n_periods

   pi = 3.1415926536_8
   twopi= 2. * pi
   omega = twopi * f_rf_Hz
   sawtooth_Cm2_vst = 0.0_8
   T_rf_s = 1./f_rf_Hz

   IF (t_s.LT.t_start_s) RETURN

   n_periods = int((t_s - t_start_s)/T_rf_s)
   phase = (t_s - t_start_s) / T_rf_s - dble(n_periods)
!  if (phase.le.asymm) then
!    arg =  pi * phase / asymm
!  else
!    arg =  pi * ( 1. + (phase-asymm) / (1.-asymm) )
!  end if

!  arg = twopi * phase**3
   arg = twopi * phase**0.5

   sawtooth_Cm2_vst =  - (j_max_Am2 / omega) * ( COS(arg) + 0.1013145 ) !zero net charge density over one period

END FUNCTION sawtooth_Cm2_vst

REAL(8) FUNCTION j_Am2_vst(t_s)
   IMPLICIT NONE
   REAL(8), PARAMETER :: j_const_Am2 = 1.9_8 !*** 1.9 A/m2
   REAL(8), PARAMETER :: t_start_s = 0.0_8
   REAL(8) t_s

   j_Am2_vst = 0.
   IF (t_s.LT.t_start_s) RETURN
   j_Am2_vst =  j_const_Am2

END FUNCTION j_Am2_vst








