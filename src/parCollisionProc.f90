!
!===================================================================================================
! Called for each colliding particle from PROCESS_COLLISIONS_WITH_ATOMS
!
SUBROUTINE COLLIDE_ELECTRON(c_kind, part_num, energy_eV)

  USE MCCollisions
  USE CurrentProblemValues
  IMPLICIT NONE

  INTEGER c_kind    ! ordering number of activated kind of collision
  INTEGER part_num  ! ordering number of particle 
  REAL(8) energy_eV ! particle energy [eV]

!print *, 'Collide_electron ', c_kind,' with energy : ',  energy_eV
  
  IF (c_kind.GT.Ncolkind_spec(1)) THEN                                  ! NULL collision, do nothing
     RETURN                                                             ! and quit the subroutine
  END IF

  SELECT CASE (Colkind_of_spec(1)%activated(c_kind))                    ! 1 = species is electron

     CASE (1)       ! Elastic, model 1

!!        CALL Add_to_stored_list(part_num)
        CALL CollideElectron_1(part_num, energy_eV) !@
        e_n_1_count = e_n_1_count + 1
        RETURN

     CASE (2)       ! Excitation, model 1

!!        CALL Add_to_stored_list(part_num)
        CALL CollideElectron_2(part_num, energy_eV) !@
        e_n_2_count = e_n_2_count + 1
        RETURN

     CASE (3)       ! Ionization, model 1

!!        CALL Add_to_stored_list(part_num)
        CALL CollideElectron_3(part_num, energy_eV) !@
        e_n_3_count = e_n_3_count + 1
        RETURN

     CASE (4)       ! Turbulence, model 1                                                            !@#$

!!        CALL Add_to_stored_list(part_num)
        CALL CollideElectron_4(part_num) !@                                                          !@#$
        e_t_4_count = e_t_4_count + 1                                                                !@#$
        RETURN                                                                                       !@#$

   END SELECT

END SUBROUTINE COLLIDE_ELECTRON

!====================================================================================================
! Called for each colliding particle from PROCESS_COLLISIONS_WITH_ATOMS
!
SUBROUTINE COLLIDE_ION(c_kind, part_num, Vx_n, Vy_n, Vz_n)

  USE MCCollisions
  USE CurrentProblemValues
  IMPLICIT NONE

  INTEGER c_kind             ! ordering number of activated kind of collision
  INTEGER part_num           ! ordering number of particle 

  REAL(8) Vx_n, Vy_n, Vz_n   ! velocity components of colliding neutral, obtained in PROCESS_COLLISIONS_WITH_ATOMS  

!print *, 'Collide_ion ', col_kind

  IF (c_kind.GT.Ncolkind_spec(2)) THEN                                  ! NULL collision, do nothing
     RETURN                                                             ! and quit the subroutine
  END IF

  SELECT CASE (Colkind_of_spec(2)%activated(c_kind))                  ! 2 = species is ions

     CASE (1)       ! Elastic, model 1

!        CALL Add_to_stored_list(part_num)
        CALL CollideIon_1(part_num, Vx_n, Vy_n, Vz_n) !@
        i_n_1_count = i_n_1_count + 1
        RETURN

     CASE (2)       ! Charge exchange, model 1

!        CALL Add_to_stored_list(part_num)
        CALL CollideIon_2(part_num, Vx_n, Vy_n, Vz_n) !@
        i_n_2_count = i_n_2_count + 1
        RETURN

     CASE (3)       ! Turbulence, model 1                                            !@#$

!        CALL Add_to_stored_list(part_num)
        CALL CollideIon_3(part_num) !@                                               !@#$
        i_t_3_count = i_t_3_count + 1                                                !@#$
        RETURN                                                                       !@#$

  END SELECT

END SUBROUTINE COLLIDE_ION

!====================================================================================================
! Frequency of electron-neutral collisions, one of all possible kinds
!
REAL(8) FUNCTION FREQUENCY_OF_COLLISION(energy_eV, active_coll, s)

  USE MCCollisions
  USE CurrentProblemValues, ONLY : delta_t_s
  IMPLICIT NONE

  REAL    energy_eV     ! particle energy [eV]
  INTEGER active_coll   ! index of collision procedure
  INTEGER s             ! species

  REAL(8) Frequency_EN_s1_1
  REAL(8) Frequency_EN_s1_2
  REAL(8) Frequency_EN_s1_3
  REAL(8) Frequency_ET_s1_4                                                           !@#$

  REAL(8) Frequency_IN_s1_1
  REAL(8) Frequency_IN_s1_2
  REAL(8) Frequency_IT_s1_3                                                           !@#$
  
  IF (s.EQ.1) THEN
     SELECT CASE (active_coll)
        CASE (1)       ! Elastic, model 1
           FREQUENCY_OF_COLLISION = Frequency_EN_s1_1(energy_eV) * delta_t_s
           RETURN
        CASE (2)       ! Excitation, model 1
           FREQUENCY_OF_COLLISION = Frequency_EN_s1_2(energy_eV) * delta_t_s
           RETURN
        CASE (3)       ! Ionization, model 1
           FREQUENCY_OF_COLLISION = Frequency_EN_s1_3(energy_eV) * delta_t_s
           RETURN
        CASE (4)       ! Turbulence, model 1                                       !@#$
           FREQUENCY_OF_COLLISION = Frequency_ET_s1_4(energy_eV) * delta_t_s       !@#$
           RETURN                                                                  !@#$
     END SELECT
  END IF

  IF (s.EQ.2) THEN
     SELECT CASE (active_coll)
        CASE (1)       ! Elastic, model 1
           FREQUENCY_OF_COLLISION = Frequency_IN_s1_1(energy_eV) * delta_t_s
           RETURN
        CASE (2)       ! Charge exchange, model 1
           FREQUENCY_OF_COLLISION = Frequency_IN_s1_2(energy_eV) * delta_t_s
           RETURN
        CASE (3)       ! Turbulence, model 1                                       !@#$
           FREQUENCY_OF_COLLISION = Frequency_IT_s1_3(energy_eV) * delta_t_s       !@#$
           RETURN                                                                  !@#$
     END SELECT
  END IF

  FREQUENCY_OF_COLLISION = 0.0_8   ! Zero if we have no implemented combination of indexes "active_coll,s" (like "2,4")

END FUNCTION FREQUENCY_OF_COLLISION 

!************************** elastic electron-neutral collisions, model 1 *****************************
! FREQUENCY  
!
REAL(8) FUNCTION Frequency_EN_s1_1(energy_eV)

  USE MCCollisions, ONLY : N_en_elast, Energy_en_elast_eV, CrSect_en_elast_m2, N_neutral_m3
  USE CurrentProblemValues, ONLY : m_e_kg, e_Cl
  IMPLICIT NONE

  REAL energy_eV, aa, bb
  REAL crsect_m2
  INTEGER j

  REAL(8) fj_s1, fjp1_s1

  IF (energy_eV.LE.1.0d-6) THEN
     Frequency_EN_s1_1 = 0.0_8
     RETURN
  END IF

  IF (energy_eV.GE.Energy_en_elast_eV(N_en_elast)) THEN

!     Frequency_EN_s1_1 = N_neutral_m3 * CrSect_en_elast_m2(N_en_elast) * SQRT(2.0_8 * energy_eV * e_Cl / m_e_kg)
!*** this should work better, 06/06/14:
     aa = energy_eV/Energy_en_elast_eV(N_en_elast)
     bb = 1./aa
     Frequency_EN_s1_1 = bb * N_neutral_m3 * CrSect_en_elast_m2(N_en_elast) * SQRT(2.0_8 * energy_eV * e_Cl / m_e_kg)
  
  ELSE 

     DO j = 1, N_en_elast - 1
       IF ((energy_eV.GE.Energy_en_elast_eV(j)).AND.(energy_eV.LT.Energy_en_elast_eV(j+1))) EXIT
     END DO
     
     j = MIN(j, N_en_elast - 1)

     IF (Energy_en_elast_eV(j).EQ.Energy_en_elast_eV(j+1)) THEN

        Frequency_EN_s1_1 = N_neutral_m3 * CrSect_en_elast_m2(j+1) * SQRT(2.0_8 * energy_eV * e_Cl / m_e_kg)
        
     ELSE
     
!!        fj_s1   = N_neutral_m3 * CrSect_en_elast_m2(j)   * SQRT(2.0_8 * Energy_en_elast_eV(j)   * e_Cl / m_e_kg)
!!        fjp1_s1 = N_neutral_m3 * CrSect_en_elast_m2(j+1) * SQRT(2.0_8 * Energy_en_elast_eV(j+1) * e_Cl / m_e_kg)

!***** LSP-like interpolation:
        fj_s1   = N_neutral_m3 * CrSect_en_elast_m2(j)   * SQRT(2.0_8 * energy_eV * e_Cl / m_e_kg)
        fjp1_s1 = N_neutral_m3 * CrSect_en_elast_m2(j+1) * SQRT(2.0_8 * energy_eV * e_Cl / m_e_kg)

!! note, we interpolate the frequency, not the cross section
! when the cross-section is interpolated, the curve frequency(energy) may become non-monotonic due to contribution of factor sqrt(energy)
        Frequency_EN_s1_1 = fj_s1 + (energy_eV - Energy_en_elast_eV(j)) * (fjp1_s1 - fj_s1) / (Energy_en_elast_eV(j+1) - Energy_en_elast_eV(j))
                                          
     END IF                                     
  
  END IF

  RETURN

END FUNCTION Frequency_EN_s1_1 
!
! COLLISION PROCESSING -------------------------------------------------------------------------------
!
SUBROUTINE CollideElectron_1(num, energy_eV)

  USE MCCollisions
  USE CurrentProblemValues
  USE ParallelOperationValues
  USE Diagnostics, ONLY : Rate_energy_coll, P_heat_cell
  USE ElectronInjection, ONLY : UseSmartTagsFlag
  USE mt19937
  IMPLICIT NONE

  INTEGER num          ! ordering number of particle 
  REAL(8) energy_eV    ! particle energy [eV]
  integer cell_number

  REAL    RAN

  REAL(8) R                ! random number             
  REAL(8) Ksi              ! scattering angle (relative to the initial direction)
  REAL(8) CosKsi, SinKsi     
  REAL(8) Fi               ! azimuthal scattering angle 
  REAL(8) CosFi, SinFi     
  REAL(8) Vx0, Vy0, Vz0
  REAL(8) Vx, Vy, Vz       ! velocity components in rest frame of the atom, before scattering
  REAL(8) Vx_s, Vy_s, Vz_s ! velocity components, after scattering
  REAL(8) V, V_xy, a, b            
  REAL(8) delta_rel         ! electron energy drop due to collision
  REAL(8) alpha            ! coefficient, accounting the electron energy drop

  REAL(8) energy_change    ! change of energy of colliding electron
  real(8) s, f_okh, Vx_n, Vy_n, Vz_n

!  IF (energy_eV.LT.0.0_8) THEN 
!     PRINT '(/2x,"Process ",i3," : Potential ERROR in CollideElectron_1 (elastic e-n collisions):")', Rank_of_process
!     PRINT  '(2x,"Particle energy ",f10.3,"(eV) is EXACT ZERO")', energy_eV
!     PRINT  '(2x,"******* ******* THIS COLLISION EVENT WILL BE SKIPPED ******* *******")'
!     e_n_1_count = e_n_1_count - 1
!     RETURN
!  END IF

! instead of a mild safety rule above, let us not scatter particles with energy below, say, 0.000001 eV
  IF (energy_eV.LT.1.0d-6) THEN
     e_n_1_count = e_n_1_count - 1
     RETURN
  END IF

   call GetMaxwellVelocity(Vx_n)
   call GetMaxwell2D(Vy_n, Vz_n)
!!   call GetMaxwellVelocity(Vy_n)
!!   call GetMaxwellVelocity(Vz_n)
   Vx_n = alpha_Vscl * Vx_n
   Vy_n = alpha_Vscl * Vy_n
   Vz_n = alpha_Vscl * Vz_n

! Calculate the scattering angle relative to the initial direction of electron velocity
  R = grnd()
! #####  CosKsi = (2.0_8 + energy_eV - 2.0_8 * (1.0_8 + energy_eV)**R) / energy_eV #####
! the formula above was in the older code and it was based on Surendra's differential cross section
! below is the corrected expression from Okhrimovsky et al., Phys.Rev.E, 65, 037402 (2002).
!  CosKsi = 1.0_8 - 2.0_8 * DBLE(R) / (1.0_8 + 8.0_8 * (energy_eV / 27.21_8) * (1.0_8 - DBLE(R)))
  SELECT CASE (Neutral_flag)
    CASE (0)  ! Helium
        s = SQRT(energy_eV)
        f_okh = 0.974 - 13.6/((s-1.82)**2 + 11.0)         ! Helium for use with cross-sections that came with EDIPIC
                                                          ! this approximates 
    CASE (1)  ! Argon
        f_okh = 0.0_8                                     ! Argon
  END SELECT
  CosKsi = dble( 1.- 2. * R * (1. - f_okh) / (1. + f_okh * (1. - 2. * R)) )

  CosKsi = MAX(MIN(0.999999999999_8, CosKsi), -0.999999999999_8)   !############ to avoid an unlikely situation when |CosKsi|>1
  Ksi = ACOS(CosKsi)
  SinKsi = SIN(Ksi)
! Calculate the azimuthal scattering angle
!  R = RANF()
  R = grnd()
  Fi = R * 6.28318530718_8
  CosFi = COS(Fi)
  SinFi = SIN(Fi)
! Take the velocity 
  Vx0 = VX_of_spec(1)%part(num)
  Vy0 = VY_of_spec(1)%part(num) 
  Vz0 = VZ_of_spec(1)%part(num)
  Vx = Vx0 - Vx_n
  Vy = Vy0 - Vy_n
  Vz = Vz0 - Vz_n
! Turn the velocity  
  V    = SQRT(Vx*Vx + Vy*Vy + Vz*Vz) 
  V_xy = SQRT(Vx*Vx + Vy*Vy)
!  IF (V_xy.GT.0.0_8) THEN                                  ! was like this, changed for a new version (below)
!     a    = SinKsi * SinFi * V / V_xy                      ! to avoid possible error when Vx, Vy are very small but non-zero
!     b    = SinKsi * CosFi * Vz / V_xy                     ! (the error may ??? appear in terms with V_x,y/V_xy)
!     Vx_s = Vx * CosKsi + Vy * a + Vx * b                  !
!     Vy_s = Vy * CosKsi - Vx * a + Vy * b                  !
!     Vz_s = Vz * CosKsi          - V_xy * SinKsi * CosFi   !
  IF (V_xy.GT.1.0d-20) THEN
     a = Vx / V_xy
     b = Vy / V_xy
     Vx_s = Vx * CosKsi + (SinFi * V * b + CosFi * Vz * a) * SinKsi
     Vy_s = Vy * CosKsi - (SinFi * V * a - CosFi * Vz * b) * SinKsi
     Vz_s = Vz * CosKsi - V_xy * CosFi * SinKsi
  ELSE
     Vx_s = ABS(Vz) * SinKsi * CosFi
     Vy_s = ABS(Vz) * SinKsi * SinFi
     Vz_s = Vz * CosKsi
  END IF

!print *, Vx_s, Vy_s, Vz_s

! Calculate the energy drop [eV]
  delta_rel = (0.001097161_8 / M_neutral_amu) * (1.0_8 - CosKsi) ! * energy_eV 
!! 0.001097161 = 2 * m_e_kg / 1_amu_kg = 2 * 9.109534e-31 / 1.660565e-27 
!! ####### NOTE, VAHEDI'S EQ.(12) EITHER HAS A MISTAKE IN THE ABOVE EXPRESSION - MISSED " * ENERGY_EV"
!! ####### OR IS GIVEN ALREADY FOR THE RELATIVE ENERGY LOSS ...

  alpha = 1.0_8 - delta_rel

! Renormalize the velocity in order to account the energy change
  IF (alpha.LT.0.0_8) THEN
     PRINT '(/2x,"Process ",i3," : ERROR in CollideElectron_1 (elastic e-n collisions):")', Rank_of_process
     PRINT  '(2x,"scattering angle cosine: CosKsi = ",e16.9)', CosKsi
     PRINT  '(2x,"factor for energy drop begore sqrt (should be positive) alpha = ",e16.9)', alpha
     PRINT  '(2x,"terminating the program")'
     STOP
  END IF
  alpha = SQRT(alpha)
  Vx_s = Vx_s * alpha + Vx_n
  Vy_s = Vy_s * alpha + Vy_n
  Vz_s = Vz_s * alpha + Vz_n
  VX_of_spec(1)%part(num) = Vx_s 
  VY_of_spec(1)%part(num) = Vy_s 
  VZ_of_spec(1)%part(num) = Vz_s

! Mark the electron as the one collided with a neutral atom ONLY if it changed the direction of its velocity
  IF (((Vx0*Vx_s).LE.0.0_8).AND.(UseSmartTagsFlag.EQ.1)) THEN
     Tag_of_spec(1)%part(num) = 0 !eTag_Coll_Neutral
  END IF

! calculate and save the change of energy
  energy_change = (Vx_s*Vx_s + Vy_s*Vy_s + Vz_s*Vz_s) - (Vx0*Vx0 + Vy0*Vy0 + Vz0*Vz0)
  Rate_energy_coll(1) = Rate_energy_coll(1) + energy_change

  cell_number = 1 + int( X_of_spec(1)%part(num) )
  if (cell_number.gt.N_cells) cell_number = N_cells
  if (cell_number.ge.1) P_heat_cell(cell_number, 1) = P_heat_cell(cell_number, 1) - energy_change

END SUBROUTINE CollideElectron_1
!=====================================================================================================

!******************* excitation (inelastic) electron-neutral collisions, model 1 *********************
! FREQUENCY 
!
REAL(8) FUNCTION Frequency_EN_s1_2(energy_eV)

  USE MCCollisions, ONLY : N_en_excit, Energy_en_excit_eV, CrSect_en_excit_m2, N_neutral_m3, Thresh_en_excit_eV
  USE CurrentProblemValues, ONLY : m_e_kg, e_Cl

  IMPLICIT NONE

  REAL energy_eV, aa, bb
  REAL crsect_m2
  INTEGER j

  REAL(8) fj_s1, fjp1_s1

  IF (energy_eV.LE.(Thresh_en_excit_eV+1.0d-6)) THEN
     Frequency_EN_s1_2 = 0.0_8
     RETURN
  END IF

  IF (energy_eV.GT.Energy_en_excit_eV(N_en_excit)) THEN

!     Frequency_EN_s1_2 = N_neutral_m3 * CrSect_en_excit_m2(N_en_excit) * SQRT(2.0_8 * energy_eV * e_Cl / m_e_kg)
! *** adjustment, 06/06/14:
     aa = energy_eV/Energy_en_excit_eV(N_en_excit)
     bb = 1./aa  
     Frequency_EN_s1_2 = bb * N_neutral_m3 * CrSect_en_excit_m2(N_en_excit) * SQRT(2.0_8 * energy_eV * e_Cl / m_e_kg)
  
  ELSE 

     DO j = 1, N_en_excit - 1
       IF ((energy_eV.GE.Energy_en_excit_eV(j)).AND.(energy_eV.LE.Energy_en_excit_eV(j+1))) EXIT
     END DO

     j = MIN(j, N_en_excit - 1)
     
     IF (Energy_en_excit_eV(j).EQ.Energy_en_excit_eV(j+1)) THEN

        Frequency_EN_s1_2 = N_neutral_m3 * CrSect_en_excit_m2(j+1) * SQRT(2.0_8 * energy_eV * e_Cl / m_e_kg)
        
     ELSE

!!        fj_s1   = N_neutral_m3 * CrSect_en_excit_m2(j)   * SQRT(2.0_8 * Energy_en_excit_eV(j)   * e_Cl / m_e_kg)
!!        fjp1_s1 = N_neutral_m3 * CrSect_en_excit_m2(j+1) * SQRT(2.0_8 * Energy_en_excit_eV(j+1) * e_Cl / m_e_kg)

!***** LSP-like interpolation instead:
        fj_s1   = N_neutral_m3 * CrSect_en_excit_m2(j)   * SQRT(2.0_8 * energy_eV * e_Cl / m_e_kg)
        fjp1_s1 = N_neutral_m3 * CrSect_en_excit_m2(j+1) * SQRT(2.0_8 * energy_eV * e_Cl / m_e_kg)

! note, we interpolate the frequency, not the cross section
! when the cross-section is interpolated, the curve frequency(energy) may become non-monotonic due to contribution of factor sqrt(energy)
        Frequency_EN_s1_2 = fj_s1 + (energy_eV - Energy_en_excit_eV(j)) * (fjp1_s1 - fj_s1) / (Energy_en_excit_eV(j+1) - Energy_en_excit_eV(j))
                                          
     END IF                                     
  
  END IF

  RETURN

END FUNCTION Frequency_EN_s1_2 
!
! COLLISION PROCESSING -------------------------------------------------------------------------------
!
SUBROUTINE CollideElectron_2(num, energy_eV) !, random_seed)

  USE ParallelOperationValues
  USE MCCollisions
  USE CurrentProblemValues
  USE Diagnostics, ONLY : Rate_energy_coll
  USE ElectronInjection, ONLY : UseSmartTagsFlag
  USE mt19937
  IMPLICIT NONE

  INTEGER num          ! ordering number of particle 
  REAL(8) energy_eV    ! initial particle energy [eV]
  REAL(8) energy_sc_eV ! particle energy after scattering [eV]

  REAL RAN

  REAL(8) R                ! random number             
  REAL(8) Ksi              ! scattering angle (relative to the initial direction)
  REAL(8) CosKsi, SinKsi     
  REAL(8) Fi               ! azimuthal scattering angle 
  REAL(8) CosFi, SinFi     
  REAL(8) Vx0, Vy0, Vz0, Vx_n, Vy_n, Vz_n
  REAL(8) Vx, Vy, Vz       ! velocity components, in neutral rest frame before scattering
  REAL(8) Vx_s, Vy_s, Vz_s ! velocity components, after scattering
  REAL(8) V, V_xy, a, b            
  REAL(8) alpha            ! coefficient, accounting the electron energy drop

  REAL(8) energy_change    ! change of energy of colliding electron
  real(8) s, f_okh, enr0_eV, eps

!  REAL(8) threshold_eV                        ! the excitation threshold, Argon, approximate
!  threshold_eV = 11.55_8                      ! the excitation threshold, Argon, approximate

! enr0_eV = 1000.

!print *, 'excitation, enter'

!  IF (energy_eV.LT.(Thresh_en_excit_eV)) THEN 
!     PRINT '(/2x,"Process ",i3," : Potential ERROR in CollideElectron_2 (excitation e-n collisions):")', Rank_of_process
!     PRINT  '(2x,"Particle energy ",f10.3,"(eV) is not more than the excitation threshold ",f6.2,"(eV)")', energy_eV, Thresh_en_excit_eV
!     PRINT  '(2x,"Such low energy particle cannot take part in this kind of collisions")'
!     PRINT  '(2x,"******* ******* THIS COLLISION EVENT WILL BE SKIPPED ******* *******")'
!     e_n_2_count = e_n_2_count - 1
!     RETURN
!!     PRINT  '(2x,"The program will be terminated now :(")'
!!     STOP
!  END IF

! instead of a mild safety rule above, let us not scatter particles with energy less than 0.000001 eV above the threshold
  IF (energy_eV.LT.(Thresh_en_excit_eV+1.0d-6)) THEN 
     e_n_2_count = e_n_2_count - 1
     RETURN
  END IF

   call GetMaxwellVelocity(Vx_n)
   call GetMaxwell2D(Vy_n, Vz_n)
!!   call GetMaxwellVelocity(Vy_n)
!!   call GetMaxwellVelocity(Vz_n)
   Vx_n = alpha_Vscl * Vx_n
   Vy_n = alpha_Vscl * Vy_n
   Vz_n = alpha_Vscl * Vz_n

! Calculate the energy of the scattered electron 
  energy_sc_eV = energy_eV - Thresh_en_excit_eV
  if (energy_sc_eV.lt.0.0_8) return

! Calculate the scattering angle relative to the initial direction of electron ! [Vahedi]: Use the modified energy "energy_sc_eV" here 
  R = grnd()
!  R = RANF()  
  
! ##### CosKsi = (2.0_8 + energy_sc_eV - 2.0_8 * (1.0_8 + energy_sc_eV)**R) / energy_sc_eV ##### 
!  CosKsi = 1.0_8 - 2.0_8 * DBLE(R) / (1.0_8 + 8.0_8 * (energy_sc_eV / 27.21_8) * (1.0_8 - DBLE(R)))

!!  s = sqrt(energy_sc_eV)
!!  f_okh = 0.974 - 13.6/((s-1.82)**2 + 11.0)
  ! f_okh = 1.

  SELECT CASE (Neutral_flag)
    CASE (0)  ! Helium
      s = SQRT(energy_eV)
      f_okh = 0.974 - 13.6/((s-1.82)**2 + 11.0)         ! Helium for use with cross-sections that came with EDIPIC
                                                        ! this approximates 
    CASE (1)  ! Argon
      f_okh = 1.0_8                                     ! Argon (should this be 1 or 0?)
  END SELECT

  CosKsi = dble( 1. - 2. * R * (1. - f_okh) / (1.+ f_okh * (1. - 2. * R)) )

  CosKsi = MAX(MIN(0.999999999999_8, CosKsi), -0.999999999999_8)   !############ to avoid an unlikely situation when |CosKsi|>1
  Ksi = ACOS(CosKsi)
  SinKsi = SIN(Ksi)
! Calculate the azimuthal scattering angle
!  R = RANF() 
  Fi = grnd() * 6.28318530718_8
  CosFi = COS(Fi)
  SinFi = SIN(Fi)
! Take the velocity 
  Vx0 = VX_of_spec(1)%part(num)
  Vy0 = VY_of_spec(1)%part(num)
  Vz0 = VZ_of_spec(1)%part(num)
  Vx = Vx0 - Vx_n
  Vy = Vy0 - Vy_n
  Vz = Vz0 - Vz_n
!print *, Vx, Vy, Vz
! Turn the velocity  
  V    = SQRT(Vx*Vx + Vy*Vy + Vz*Vz) 
  V_xy = SQRT(Vx*Vx + Vy*Vy)
!  IF (V_xy.GT.0.0_8) THEN
!     a    = SinKsi * SinFi * V / V_xy
!     b    = SinKsi * CosFi * Vz / V_xy 
!     Vx_s = Vx * CosKsi + Vy * a + Vx * b
!     Vy_s = Vy * CosKsi - Vx * a + Vy * b
!     Vz_s = Vz * CosKsi          - V_xy * SinKsi * CosFi
  IF (V_xy.GT.1.0d-20) THEN
     a = Vx / V_xy
     b = Vy / V_xy
     Vx_s = Vx * CosKsi + (SinFi * V * b + CosFi * Vz * a) * SinKsi
     Vy_s = Vy * CosKsi - (SinFi * V * a - CosFi * Vz * b) * SinKsi
     Vz_s = Vz * CosKsi - V_xy * CosFi * SinKsi
  ELSE
     Vx_s = ABS(Vz) * SinKsi * CosFi
     Vy_s = ABS(Vz) * SinKsi * SinFi
     Vz_s = Vz * CosKsi
  END IF

! Calculate the energy drop 
!  alpha    = SQRT(1.0_8 - Thresh_en_excit_eV / energy_eV)
  alpha    = SQRT(energy_sc_eV / energy_eV)
! Renormalize the velocity in order to account the energy drop
  Vx_s = Vx_s * alpha + Vx_n   
  Vy_s = Vy_s * alpha + Vy_n   
  Vz_s = Vz_s * alpha + Vz_n     
  VX_of_spec(1)%part(num) = Vx_s 
  VY_of_spec(1)%part(num) = Vy_s 
  VZ_of_spec(1)%part(num) = Vz_s 

! Mark the electron as the one collided with a neutral atom ONLY if it changed the direction of its velocity 
  IF (((Vx0*Vx_s).LE.0.0_8).AND.(UseSmartTagsFlag.EQ.1)) THEN
     Tag_of_spec(1)%part(num) = 0 !eTag_Coll_Neutral
  END IF

! calculate and save the change of energy
  energy_change = (Vx_s*Vx_s + Vy_s*Vy_s + Vz_s*Vz_s) - (Vx0*Vx0 + Vy0*Vy0 + Vz0*Vz0)
  Rate_energy_coll(1) = Rate_energy_coll(1) + energy_change

!print *, 'excitation, exit'

END SUBROUTINE CollideElectron_2
!=====================================================================================================

!******************** ionization (inelastic) electron-neutral collisions, model 1 ********************
! FREQUENCY 
!
REAL(8) FUNCTION Frequency_EN_s1_3(energy_eV)

  USE MCCollisions, ONLY : N_en_ioniz, Energy_en_ioniz_eV, CrSect_en_ioniz_m2, N_neutral_m3, Thresh_en_ioniz_eV
  USE CurrentProblemValues, ONLY : m_e_kg, e_Cl

  IMPLICIT NONE

  REAL energy_eV, aa, bb
  REAL crsect_m2
  INTEGER j

  REAL(8) fj_s1, fjp1_s1

  IF (energy_eV.LE.(Thresh_en_ioniz_eV+1.0d-6)) THEN
     Frequency_EN_s1_3 = 0.0_8
     RETURN
  END IF

  IF (energy_eV.GT.Energy_en_ioniz_eV(N_en_ioniz)) THEN

!     Frequency_EN_s1_3 = N_neutral_m3 * CrSect_en_ioniz_m2(N_en_ioniz) * SQRT(2.0_8 * energy_eV * e_Cl / m_e_kg)
! adjustment 06/06/14:
     aa = energy_eV/Energy_en_ioniz_eV(N_en_ioniz)
     bb = 1./aa
     Frequency_EN_s1_3 = bb * N_neutral_m3 * CrSect_en_ioniz_m2(N_en_ioniz) * SQRT(2.0_8 * energy_eV * e_Cl / m_e_kg)
  
  ELSE 

     DO j = 1, N_en_ioniz - 1
       IF ((energy_eV.GE.Energy_en_ioniz_eV(j)).AND.(energy_eV.LE.Energy_en_ioniz_eV(j+1))) EXIT
     END DO

     j = MIN(j, N_en_ioniz - 1)
     
     IF (Energy_en_ioniz_eV(j).EQ.Energy_en_ioniz_eV(j+1)) THEN
     
        Frequency_EN_s1_3 = N_neutral_m3 * CrSect_en_ioniz_m2(j+1) * SQRT(2.0_8 * energy_eV * e_Cl / m_e_kg)
        
     ELSE
          
!!        fj_s1   = N_neutral_m3 * CrSect_en_ioniz_m2(j)   * SQRT(2.0_8 * Energy_en_ioniz_eV(j)   * e_Cl / m_e_kg)
!!        fjp1_s1 = N_neutral_m3 * CrSect_en_ioniz_m2(j+1) * SQRT(2.0_8 * Energy_en_ioniz_eV(j+1) * e_Cl / m_e_kg)

!***** LSP-like interpolation instead:
        fj_s1   = N_neutral_m3 * CrSect_en_ioniz_m2(j)   * SQRT(2.0_8 * energy_eV * e_Cl / m_e_kg)
        fjp1_s1 = N_neutral_m3 * CrSect_en_ioniz_m2(j+1) * SQRT(2.0_8 * energy_eV * e_Cl / m_e_kg)

! note, we interpolate the frequency, not the cross section
! when the cross-section is interpolated, the curve frequency(energy) may become non-monotonic due to contribution of factor sqrt(energy)
        Frequency_EN_s1_3 = fj_s1 + (energy_eV - Energy_en_ioniz_eV(j)) * (fjp1_s1 - fj_s1) / (Energy_en_ioniz_eV(j+1) - Energy_en_ioniz_eV(j))
                                           
     END IF                                     
  
  END IF

  RETURN

END FUNCTION Frequency_EN_s1_3 
!
! COLLISION PROCESSING -------------------------------------------------------------------------------
!
SUBROUTINE CollideElectron_3(num, energy_inc_eV) 

  USE ParallelOperationValues
  USE MCCollisions
  USE CurrentProblemValues
  USE Diagnostics, ONLY : Rate_energy_coll, N_new_cell
  USE ElectronInjection, ONLY : UseSmartTagsFlag
  USE mt19937
  IMPLICIT NONE

  INTEGER num           ! ordering number of particle 
  REAL(8) energy_inc_eV ! incident electron energy [eV]

  REAL RAN

  INTEGER ALLOC_ERR

  REAL(8) energy_sc_eV     ! energy of incident electron after scattering [eV]
  REAL(8) Ksi_sc
  REAL(8) CosKsi_sc, SinKsi_sc     
  REAL(8) Fi_sc
  REAL(8) CosFi_sc, SinFi_sc     

  REAL(8) energy_ej_eV     ! energy of ejected electron [eV] 
  REAL(8) Ksi_ej
  REAL(8) CosKsi_ej, SinKsi_ej     
  REAL(8) Fi_ej
  REAL(8) CosFi_ej, SinFi_ej     

  REAL(8) V, V_xy, a, b            

  REAL(8) Vx, Vy, Vz          ! velocity components of incident electron, before scattering
  REAL(8) Vx_sc, Vy_sc, Vz_sc ! velocity components of incident electron, after scattering
  REAL(8) Vx_ej, Vy_ej, Vz_ej ! velocity components of ejected electron
  REAL(8) Vx_i, Vy_i, Vz_i    ! velocity components of produced ion

  REAL(8) alpha            

  REAL(8) R                  ! random number             

  REAL(8) energy_change_e    ! change of kinetic energy of electrons, incident and ejected
  REAL(8) energy_change_i    ! change of kinetic energy due to produced ion

  REAL B_of_Einc_eV                       ! some known function of incident electron energy
  real(8) s, f_okh, eps, enr0_eV !screening parameter, can be adjusted because diff. scattering for ionization is not quite known

  INTEGER left_node, right_node
  real(8) x0, xleft, xright, Ex0

!  REAL threshold_eV                       ! the ionization threshold

  x0 = X_of_spec(1)%part(num)

  if (x0.lt.0.0_8.or.x0.gt.dble(N_cells)) return !safety pre-caution

  SELECT CASE (Neutral_flag)
    CASE (0)  ! Helium
        B_of_Einc_eV = 15.8_8 !** Helium (Opal, Peterson and Beatty)
    CASE (1)  ! Argon
        B_of_Einc_eV = 10.0_8 ! is constant for Einc < 70 eV, Argon, approximate
  END SELECT

!  enr0_eV = 150.0_8 !screening parameter value utilized by Johan in ionizing collisions

!############################################################################################      ##    ##
!###################### NOTE: MUST BE: "B_of_Einc_eV = 8.7" ###################################      ##    ##
!############################################################################################      ##    ##
!  threshold_eV = 15.0_8                   ! the ionization threshold, Argon, approximate

!print *, 'ionization'

!  IF (energy_inc_eV.LE.Thresh_en_ioniz_eV) THEN 
!     PRINT '(/2x,"Process ",i3," : Potential ERROR in CollideElectron_3 (ionization e-n collision):")', Rank_of_process
!     PRINT  '(2x,"Particle energy ",f10.3,"(eV) is not more than the ionization threshold ",f6.2,"(eV)")',&
!                                                                    & energy_inc_eV, Thresh_en_ioniz_eV
!     PRINT  '(2x,"Such low energy particle cannot take part in this kind of collisions")'
!     PRINT  '(2x,"******* ******* THIS COLLISION EVENT WILL BE SKIPPED ******* *******")'
!     e_n_3_count = e_n_3_count - 1
!     RETURN
!!     PRINT  '(2x,"The program will be terminated now :(")'
!!     STOP
!  END IF

  IF (energy_inc_eV.LT.(Thresh_en_ioniz_eV+1.0d-6)) THEN 
     e_n_3_count = e_n_3_count - 1
     RETURN
  END IF

! Calculate the energy of the ejected electron
  R = grnd()
!  R = RANF()
  energy_ej_eV = B_of_Einc_eV * TAN(R * ATAN( 0.5_8 * (energy_inc_eV - Thresh_en_ioniz_eV) / B_of_Einc_eV ))  
! note, theoretically, for R<=1, the value of energy_ej_eV should never exceed (energy_inc_eV - Thresh_en_ioniz_eV)/2
! but below we shall enforce this manually

! I added the line below to be sure that the energy of ejected electron is always
! (a) non-zero and (b) below half of energy difference [incident - threshold]
! note, 5d-7 = 1d-6/2, where 1d-6 is the minimal positive energy difference [incident - threshold]
  energy_ej_eV = MAX(5.0d-7, MIN(energy_ej_eV, 0.5_8 * (energy_inc_eV - Thresh_en_ioniz_eV)))

!  if (energy_ej_eV.gt.(0.5_8 * (energy_inc_eV - Thresh_en_ioniz_eV))) then
!     PRINT '(/2x,"Process ",i3," : Error in CollideElectron_3:")', Rank_of_process
!     PRINT  '(2x,"energy of  ",f10.3,"(eV) is greater than ",f6.2,"(eV)")', energy_ej_eV, &
!                                              & 0.5_8 * (energy_inc_eV - Thresh_en_ioniz_eV)
!     PRINT  '(2x,"The energy conservation law can be violated")'
!     PRINT  '(2x,"The program will be terminated now :(")'
!     STOP
!  END IF

! Calculate the energy of the scattered (incident) electron
  energy_sc_eV = energy_inc_eV - Thresh_en_ioniz_eV - energy_ej_eV

!##### it is expected that here neither energy_ej_eV nor energy_sc_eV are zeros (both should be not less than 5d-7 eV) #####

!print *, energy_inc_eV, energy_sc_eV, energy_ej_eV

! Calculate the scattering angle Ksi for the incident electron ! [Vahedi]: use the modified energy "energy_sc_eV" here
 R = grnd()
!####  CosKsi_sc = (2.0_8 + energy_sc_eV - 2.0_8 * (1.0_8 + energy_sc_eV)**R) / energy_sc_eV ####
!  CosKsi_sc = 1.0_8 - 2.0_8 * DBLE(R) / (1.0_8 + 8.0_8 * (energy_sc_eV / 27.21_8) * (1.0_8 - DBLE(R))) !new default from the author

!! s = sqrt(energy_sc_eV)
!! f_okh = 0.974 - 13.6/((s-1.82)**2 + 11.0)
  ! f_okh = 1.
! eps = 4. * energy_sc_eV / enr0_eV !screened Coulomb, adjustable
! f_okh = eps / (1. + eps) 

  SELECT CASE (Neutral_flag)
    CASE (0)  ! Helium
      s = SQRT(energy_inc_eV)
      f_okh = 0.974 - 13.6/((s-1.82)**2 + 11.0)         ! Helium for use with cross-sections that came with EDIPIC
                                                        ! this approximates 
    CASE (1)  ! Argon
      f_okh = 1.0_8                                     ! Argon (should this be 1 or 0?)
  END SELECT

 CosKsi_sc = dble( 1.- 2. * R * (1. - f_okh)/(1. + f_okh * (1. - 2. * R)) )

! CosKsi_sc = sqrt( energy_sc_eV / (energy_inc_eV - Thresh_en_ioniz_eV) ) ! *** a la Boeuf ***

  CosKsi_sc = MAX(MIN(0.999999999999_8, CosKsi_sc), -0.999999999999_8)   !############ to avoid an unlikely situation when |CosKsi_sc|>1
  Ksi_sc    = ACOS(CosKsi_sc)
  SinKsi_sc = SIN(Ksi_sc)

!  SinKsi_sc = sqrt (1. - CosKsi_sc**2) ! a la Boeuf

! Calculate the azimuthal scattering angle for the incident electron
  R = grnd()
  Fi_sc    = R * 6.28318530718_8
  CosFi_sc = COS(Fi_sc)
  SinFi_sc = SIN(Fi_sc)

! Calculate the scattering angle Ksi for the ejected electron
 R = grnd()
!####  CosKsi_ej = (2.0_8 + energy_ej_eV - 2.0_8 * (1.0_8 + energy_ej_eV)**R) / energy_ej_eV ####
!  CosKsi_ej = 1.0_8 - 2.0_8 * DBLE(R) / (1.0_8 + 8.0_8 * (energy_ej_eV / 27.21_8) * (1.0_8 - DBLE(R)))

!!  s = sqrt(energy_ej_eV)
!!  f_okh = 0.974 - 13.6/((s-1.82)**2 + 11.0)
  f_okh = 0.0_8
!  eps = 4. * energy_ej_eV / enr0_eV
!  f_okh = eps / (1. + eps)

  CosKsi_ej = dble( 1.- 2.* R * (1. - f_okh) / (1.+ f_okh * (1. - 2.* R)) )
  CosKsi_ej = MAX(MIN(0.999999999999_8, CosKsi_ej), -0.999999999999_8)   !############ to avoid an unlikely situation when |CosKsi_ej|>1
  Ksi_ej    = ACOS(CosKsi_ej)
  SinKsi_ej = SIN(Ksi_ej)
  
!  CosKsi_ej =  SinKsi_sc ! a la Boeuf, co-planar scattering
!  SinKsi_ej = -CosKsi_sc

! Calculate the azimuthal scattering angle for the ejected electron
!  R = RANF()
  Fi_ej    = grnd() * 6.28318530718_8
  CosFi_ej = COS(Fi_ej)
  SinFi_ej = SIN(Fi_ej)
!    CosFi_ej = CosFi_sc ! a la Boeuf, co-planar scattering
!    SinFi_ej = SinFi_sc

! Take the velocity of the incident electron before the scattering
  Vx = VX_of_spec(1)%part(num)
  Vy = VY_of_spec(1)%part(num)
  Vz = VZ_of_spec(1)%part(num)

! Precalculate   
  V    = SQRT(Vx*Vx + Vy*Vy + Vz*Vz) 
  V_xy = SQRT(Vx*Vx + Vy*Vy)

! Calculate the velocity of the INCIDENT ELECTRON after scattering (turn the velocity)
!  IF (V_xy.GT.0.0_8) THEN
!     a     = SinKsi_sc * SinFi_sc * V / V_xy
!     b     = SinKsi_sc * CosFi_sc * Vz / V_xy 
!     Vx_sc = Vx * CosKsi_sc + Vy * a + Vx * b
!     Vy_sc = Vy * CosKsi_sc - Vx * a + Vy * b
!     Vz_sc = Vz * CosKsi_sc - V_xy * SinKsi_sc * CosFi_sc
  IF (V_xy.GT.1.0d-20) THEN
     a = Vx / V_xy
     b = Vy / V_xy
     Vx_sc = Vx * CosKsi_sc + (SinFi_sc * V * b + CosFi_sc * Vz * a) * SinKsi_sc
     Vy_sc = Vy * CosKsi_sc - (SinFi_sc * V * a - CosFi_sc * Vz * b) * SinKsi_sc
     Vz_sc = Vz * CosKsi_sc - V_xy * CosFi_sc * SinKsi_sc
  ELSE
     Vx_sc = ABS(Vz) * SinKsi_sc * CosFi_sc
     Vy_sc = ABS(Vz) * SinKsi_sc * SinFi_sc
     Vz_sc = Vz * CosKsi_sc
  END IF

! Calculate the relative energy drop 
  alpha = SQRT(energy_sc_eV / energy_inc_eV)
! Renormalize the velocity of incident electron in order to account the energy drop
  Vx_sc = Vx_sc * alpha
  Vy_sc = Vy_sc * alpha
  Vz_sc = Vz_sc * alpha
  VX_of_spec(1)%part(num) = Vx_sc
  VY_of_spec(1)%part(num) = Vy_sc
  VZ_of_spec(1)%part(num) = Vz_sc

! Mark the incident electron as the one collided with a neutral atom ONLY if it changed the direction of its velocity
  IF (((Vx*Vx_sc).LE.0.0_8).AND.(UseSmartTagsFlag.EQ.1)) THEN
     Tag_of_spec(1)%part(num) = 0 !eTag_Coll_Neutral
  END IF

! Calculate the velocity of the EJECTED ELECTRON after scattering (turn the velocity)
!  IF (V_xy.GT.0.0_8) THEN
!     a     = SinKsi_ej * SinFi_ej * V / V_xy
!     b     = SinKsi_ej * CosFi_ej * Vz / V_xy 
!     Vx_ej = Vx * CosKsi_ej + Vy * a + Vx * b
!     Vy_ej = Vy * CosKsi_ej - Vx * a + Vy * b
!     Vz_ej = Vz * CosKsi_ej          - V_xy * SinKsi_ej * CosFi_ej
  IF (V_xy.GT.1.0d-20) THEN
     a = Vx / V_xy
     b = Vy / V_xy
     Vx_ej = Vx * CosKsi_ej + (SinFi_ej * V * b + CosFi_ej * Vz * a) * SinKsi_ej
     Vy_ej = Vy * CosKsi_ej - (SinFi_ej * V * a - CosFi_ej * Vz * b) * SinKsi_ej
     Vz_ej = Vz * CosKsi_ej - V_xy * CosFi_ej * SinKsi_ej
  ELSE
     Vx_ej = ABS(Vz) * SinKsi_ej * CosFi_ej
     Vy_ej = ABS(Vz) * SinKsi_ej * SinFi_ej
     Vz_ej = Vz * CosKsi_ej
  END IF

!print *, Vx**2 + Vy**2 + Vz**2, Vx_ej**2 + Vy_ej**2 + Vz_ej**2   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calculate the relative energy drop 
  alpha = SQRT(energy_ej_eV / energy_inc_eV)
! Renormalize the velocity of EJECTED ELECTRON in order to account the energy drop
  Vx_ej = Vx_ej * alpha
  Vy_ej = Vy_ej * alpha
  Vz_ej = Vz_ej * alpha

!print *, energy_inc_eV, ((Vx_sc**2+Vy_sc**2+Vz_sc**2)*T_e_eV/N_box_vel**2)/energy_sc_eV, ((Vx_ej**2+Vy_ej**2+Vz_ej**2)*T_e_eV/N_box_vel**2)/energy_ej_eV
     left_node = INT(x0)
     right_node = left_node + 1
     xleft = dble(left_node)
     xright = dble(right_node)
     if (right_node.GT.dble(N_cells)) right_node = N_cells
     Ex0 = EX(right_node) * (x0 -xleft) + EX(left_node) * (xright - x0)
     
! Save the parameters of the ejected electron in the linked list
  IF (ASSOCIATED(Current_electron)) THEN
     N_inject(1) = N_inject(1) + 1
!     Current_electron%X       = X_of_spec(1)%part(num)
     Current_electron%X       = x0
     Current_electron%VX      = Vx_ej
     Current_electron%VY      = Vy_ej
     Current_electron%VZ      = Vz_ej
     Current_electron%AX      = 0.0_8
!     Current_electron%AX      = QMs(1) * Ex0
     Current_electron%Tag     = eTag_Coll_Neutral                ! mark the ejected electron as the one collided with a neutral atom
     ALLOCATE(Current_electron%next, STAT = ALLOC_ERR)
     IF (ALLOC_ERR.NE.0) THEN
        PRINT '(/2x,"Process ",i3," : Error in ALLOCATE Current_electron%next !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
     Current_electron => Current_electron%next
     NULLIFY(Current_electron%next)
  ELSE
     PRINT '(/2x,"Process ",i3," : Error in CollideElectron_3:")', Rank_of_process
     PRINT  '(2x,"Current_electron is NOT associated!")'
     PRINT  '(2x,"The program will be terminated now :(")'
     STOP
  END IF

! calculate and save the change of energy for electrons
  energy_change_e = (Vx_ej**2 + Vy_ej**2 + Vz_ej**2) + &    ! energy of ejected electron PLUS
                  & (Vx_sc**2 + Vy_sc**2 + Vz_sc**2) - &    ! energy of incident electron AFTER collision (scattered) MINUS
                  & (Vx**2    + Vy**2    + Vz**2)           ! energy of incident electron BEFORE collision
  Rate_energy_coll(1) = Rate_energy_coll(1) + energy_change_e

 ! Take random velocity from normalized maxwell distribution  
  CALL GetMaxwellVelocity(Vx_i) 
  CALL GetMaxwell2D(Vy_i, Vz_i)
!!  CALL GetMaxwellVelocity(Vy_i) 
!!  CALL GetMaxwellVelocity(Vz_i)

!    CALL GetMaxwellVector( Vx_i, Vy_i, Vz_i ) 

! Use the factor above to obtain the dim-less velocity (V * N_box_vel / V_te) of the produced ion
  Vx_i = Vx_i * alpha_Vscl
  Vy_i = Vy_i * alpha_Vscl
  Vz_i = Vz_i * alpha_Vscl

! save the parameters of the produced ion in the linked list
  IF (ASSOCIATED(Current_ion)) THEN
     N_inject(2) = N_inject(2) + 1
!     Current_ion%X       = X_of_spec(1)%part(num)
     Current_ion%X       = x0
     Current_ion%VX      = Vx_i
     Current_ion%VY      = Vy_i
     Current_ion%VZ      = Vz_i
     Current_ion%AX      = 0.0_8
!     Current_ion%AX      = QMs(2) * Ex0
     Current_ion%Tag     = 0                            ! this value can be modified here later !!!!!!!!!!!
     ALLOCATE(Current_ion%next, STAT = ALLOC_ERR)
     IF (ALLOC_ERR.NE.0) THEN
        PRINT '(/2x,"Process ",i3," : Error in ALLOCATE Current_ion%next !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
     Current_ion => Current_ion%next
     NULLIFY(Current_ion%next)
  ELSE
     PRINT '(/2x,"Process ",i3," : Error in CollideElectron_3:")', Rank_of_process
     PRINT  '(2x,"Current_ion is NOT associated!")'
     PRINT  '(2x,"The program will be terminated now :(")'
     STOP
  END IF

! calculate and save the change of energy for ion
  energy_change_i = Vx_i*Vx_i + Vy_i*Vy_i + Vz_i*Vz_i        
  Rate_energy_coll(2) = Rate_energy_coll(2) + energy_change_i

! account for the density of emitted particle (must be performed in order to obtain correct particles density in diagnostics)
  left_node  = INT(X_of_spec(1)%part(num))
  right_node = left_node + 1
  Q_strm_spec(left_node, 1)  = Q_strm_spec(left_node, 1)  + right_node - X_of_spec(1)%part(num)    ! electrons
  Q_strm_spec(right_node, 1) = Q_strm_spec(right_node, 1) + X_of_spec(1)%part(num) - left_node     ! 
  Q_strm_spec(left_node, 2)  = Q_strm_spec(left_node, 2)  + right_node - X_of_spec(1)%part(num)    ! ions
  Q_strm_spec(right_node, 2) = Q_strm_spec(right_node, 2) + X_of_spec(1)%part(num) - left_node     !

!******* counter of newly created i-e pairs for diagnostics, to be cleared each time after writing into the file:
  N_new_cell(right_node) = N_new_cell(right_node) + 1
  
END SUBROUTINE CollideElectron_3
!=====================================================================================================

!**************************** electron turbulence collisions, model 1 ********************************
! FREQUENCY  
!
REAL(8) FUNCTION Frequency_ET_s1_4(energy_eV)                 !@#$

  USE MCCollisions
  IMPLICIT NONE

  REAL energy_eV

  Frequency_ET_s1_4 = Freq_turb_e_1_s1

  RETURN

END FUNCTION Frequency_ET_s1_4 
!
! COLLISION PROCESSING -------------------------------------------------------------------------------
!
SUBROUTINE CollideElectron_4(num)       !@#$

  USE MCCollisions
  USE CurrentProblemValues          
  USE mt19937
  IMPLICIT NONE

  INTEGER num          ! ordering number of particle 

  REAL    RAN

  REAL(8) R                ! random number             
  REAL(8) Fi               ! azimuthal scattering angle (in y-z-plane)
  REAL(8) CosFi, SinFi     
  REAL(8) Vy, Vz           ! velocity components, before scattering

! Calculate the azimuthal scattering angle
  R = grnd()
  Fi = R * 6.28318530718_8
  CosFi = COS(Fi)
  SinFi = SIN(Fi)

! Take the velocity components, perpendicular to the external X-magnetic field
  Vy = VY_of_spec(1)%part(num)
  Vz = VZ_of_spec(1)%part(num)

! Assign the scattered (turned) perpendicular velocity components
  VY_of_spec(1)%part(num) = Vy * CosFi - Vz * SinFi
  VZ_of_spec(1)%part(num) = Vy * SinFi + Vz * CosFi

! Later (i.e. in further diagnostics) the value of the tag can be modified here to mark the collided particle

! energy is not changed here

END SUBROUTINE CollideElectron_4
!=====================================================================================================

!****************************** elastic ion-neutral collisions, model 1 ******************************
! FREQUENCY
!
REAL(8) FUNCTION Frequency_IN_s1_1(energy_eV)

  USE MCCollisions
  IMPLICIT NONE

  REAL energy_eV

  Frequency_IN_s1_1 = 1.0e7
  RETURN

END FUNCTION Frequency_IN_s1_1 
!
! COLLISION PROCESSING -------------------------------------------------------------------------------
!
SUBROUTINE CollideIon_1(num, Vx_n, Vy_n, Vz_n) !, random_seed)
  USE MCCollisions
  USE CurrentProblemValues
  USE Diagnostics, ONLY : Rate_energy_coll
  USE mt19937
  IMPLICIT NONE

  INTEGER num                 ! ordering number of particle 
  REAL(8) Vx_n, Vy_n, Vz_n    ! velocity components of collided neutral atom

  REAL RAN
  REAL R                   ! random number             

  REAL(8) Ksi              ! scattering angle (relative to the initial direction)
  REAL(8) CosKsi, SinKsi     
  REAL(8) Fi               ! azimuthal scattering angle 
  REAL(8) CosFi, SinFi     
  REAL(8) Vx, Vy, Vz       ! velocity components, before scattering
  REAL(8) Vx_s, Vy_s, Vz_s ! velocity components, after scattering
  REAL(8) V, V_xy, a, b, Vx_cm, Vy_cm, Vz_cm           

  REAL(8) energy_change    ! change of energy of colliding ion

! We assume that the masses of ion and neutral atom are equal, 
! In the C.O.M. frame scattering is isotropic.

! Calculate the scattering angle Ksi for the ion
  CosKsi = 2.0_8 * grnd() - 1.0_8
  SinKsi = sqrt((1.0_8 - CosKsi) * (1.0_8 + CosKsi))

! Calculate the azimuthal scattering angle for the ion
  Fi    = grnd() * 6.28318530718_8
  CosFi = COS(Fi)
  SinFi = SIN(Fi)

! Take the initial ion velocity and transform:
  Vx   = VX_of_spec(2)%part(num) !- Vx_n
  Vy   = VY_of_spec(2)%part(num) !- Vy_n
  Vz   = VZ_of_spec(2)%part(num) !- Vz_n
  energy_change = Vx**2 + Vy**2 + Vz**2

  Vx_cm = 0.5 * (Vx + Vx_n)
  Vy_cm = 0.5 * (Vy + Vy_n)
  Vz_cm = 0.5 * (Vz + Vz_n)

  Vx = Vx - Vx_cm
  Vy = Vy - Vy_cm
  Vz = Vz - Vz_cm
  
! Scattering the ion in the CM frame.
! Turn the velocity in that frame, using Ksi and Fi. After that, transform the velocity to the laboratory frame.
  
  V = SQRT(Vx * Vx + Vy * Vy + Vz * Vz)  
  Vx_s = V * CosKsi
  Vy_s = V * SinKsi * SinFi
  Vz_s = V * SinKsi * CosFi
!  Vx_n = -Vx_s + Vx_cm
!  Vy_n = -Vy_s + Vy_cm
!  Vz_n = -Vz_s + Vz_cm
  VX_of_spec(2)%part(num) = Vx_s + Vx_cm
  VY_of_spec(2)%part(num) = Vy_s + Vy_cm
  VZ_of_spec(2)%part(num) = Vz_s + Vz_cm
  
! the tag of collided particle can be modified here !!!!!!!!!!!!!!!
  
! calculate and save the change of energy
  energy_change = VX_of_spec(2)%part(num)**2 + VY_of_spec(2)%part(num)**2 + VZ_of_spec(2)%part(num)**2 - energy_change
  Rate_energy_coll(2) = Rate_energy_coll(2) + energy_change

 END SUBROUTINE CollideIon_1
!=====================================================================================================

!************************* charge exchange ion-neutral collisions, model 1 ***************************
! FREQUENCY 
!
!REAL(8) FUNCTION Frequency_IN_s1_2(energy_eV)

!  USE MCCollisions
!  IMPLICIT NONE

!  REAL energy_eV

!  Frequency_IN_s1_2 = 1.0e7

!END FUNCTION Frequency_IN_s1_2 

!******************************
REAL(8) FUNCTION Frequency_IN_s1_2(energy_eV)
! energy corresponds to relative velocity
  USE MCCollisions
  USE CurrentProblemValues, ONLY : e_Cl, amu_kg
  IMPLICIT NONE
  REAL energy_eV, V_rel_ms, sigma_ct_m2, sigma_ct_m2_1eV, gammaR, factor, energy_eff_eV, estar_eV, sigma_star_m2
  gammaR = 10.

  SELECT CASE (Neutral_flag)
    CASE (0)  ! Helium
      sigma_ct_m2_1eV = 2.8e-19 !!** Helium gas
    CASE (1)  ! Argon
      sigma_ct_m2_1eV = 5.53e-19 !!** Argon gas
  END SELECT
  
! charge exchange cross-section for helium, per BM Smirnov:
    energy_eff_eV = energy_eV
    if (energy_eV .le. 1.e-4) energy_eff_eV = 1.e-4
    SELECT CASE (Neutral_flag)
      CASE (0)  ! Helium
        factor = ( 1. - (1./gammaR) * 0.5 * log(energy_eff_eV) ) ** 2 !!** Helium gas
      CASE (1)  ! Argon
        factor = (1.0_8 - 0.0543 * log(energy_eff_eV)) ** 2 !!** Argon gas, per Maiorov
    END SELECT
    sigma_ct_m2 = sigma_ct_m2_1eV * factor
    if (sigma_ct_m2 .ge. 1.e-18) sigma_ct_m2 = 1.e-18
    V_rel_ms = SQRT( (2. * e_Cl * energy_eV ) / (M_neutral_amu * amu_kg ) )
!*!  sigma_ct_m2 = 3.E-19 !*!
   Frequency_IN_s1_2 = dble(V_rel_ms * N_neutral_m3 * sigma_ct_m2)

END FUNCTION Frequency_IN_s1_2
!*****************************
!
! COLLISION PROCESSING -------------------------------------------------------------------------------
!
SUBROUTINE CollideIon_2( num, Vx_n, Vy_n, Vz_n )

  USE MCCollisions
  USE CurrentProblemValues
  USE Diagnostics, ONLY : Rate_energy_coll, P_heat_cell
  IMPLICIT NONE

  INTEGER num               ! ordering number of particle 
  REAL(8) Vx_n, Vy_n, Vz_n  ! velocity components of neutral atom
  integer cell_number

  REAL(8) energy_change    ! change of energy of colliding electron

! Now we do not follow the neutral AFTER collision ...  

! remember the initial ion energy
  energy_change = VX_of_spec(2)%part(num)**2 + VY_of_spec(2)%part(num)**2 + VZ_of_spec(2)%part(num)**2 

! Set the ion velocity equal to the velocity of atom (we assume that neutral and ion just exchange one electron)
  VX_of_spec(2)%part(num) = Vx_n
  VY_of_spec(2)%part(num) = Vy_n
  VZ_of_spec(2)%part(num) = Vz_n

! the tag of collided particle can be modified here  !!!!!!!!!!!!!!

! calculate and save the change of energy
  energy_change = Vx_n**2 + Vy_n**2 + Vz_n**2 - energy_change
  Rate_energy_coll(2) = Rate_energy_coll(2) + energy_change
  cell_number = 1 + int( X_of_spec(2)%part(num) )
  if (cell_number.gt.N_cells) cell_number = N_cells
  if (cell_number.ge.1) P_heat_cell(cell_number, 2) = P_heat_cell(cell_number, 2) - energy_change

END SUBROUTINE CollideIon_2
!=====================================================================================================

!******************************* ion turbulence collisions, model 1 **********************************
! FREQUENCY  
!
REAL(8) FUNCTION Frequency_IT_s1_3(energy_eV)                 !@#$

  USE MCCollisions
  IMPLICIT NONE

  REAL energy_eV

  Frequency_IT_s1_3 = Freq_turb_i_1_s1

  RETURN

END FUNCTION Frequency_IT_s1_3 
!
! COLLISION PROCESSING -------------------------------------------------------------------------------
!
SUBROUTINE CollideIon_3(num)       !@#$

  USE MCCollisions
  USE CurrentProblemValues          
  USE mt19937
  IMPLICIT NONE

  INTEGER num          ! ordering number of particle 

  REAL    RAN

!  INTEGER indx             ! current block index
  REAL    R                ! random number             
  REAL(8) Fi               ! azimuthal scattering angle (in y-z-plane)
  REAL(8) CosFi, SinFi     
  REAL(8) Vy, Vz           ! velocity components, before scattering

! Calculate the azimuthal scattering angle
  R = grnd()
  Fi = R * 6.28318530718_8
  CosFi = COS(Fi)
  SinFi = SIN(Fi)

! Take the velocity components, perpendicular to the external X-magnetic field
  Vy = VY_of_spec(2)%part(num)
  Vz = VZ_of_spec(2)%part(num)

! Assign the scattered (turned) perpendicular velocity components
  VY_of_spec(2)%part(num) = Vy * CosFi - Vz * SinFi
  VZ_of_spec(2)%part(num) = Vy * SinFi + Vz * CosFi

! the tag of collided particle can be modified here !!!!!!!!!!!!!!!

! energy is not changed here

END SUBROUTINE CollideIon_3
!=====================================================================================================
