SUBROUTINE SAVE_CHECKPOINT

   USE CurrentProblemValues
   USE ParallelOperationValues
   USE Diagnostics
   USE Snapshots
   USE mt19937

   IMPLICIT NONE

   INTEGER n ! number of particle

   REAL R

   IF (T_cntr.NE.Save_check_Tcntr) RETURN

   IF (Rank_of_process.EQ.0) THEN

      OPEN (9, FILE = check_g_filename, STATUS = 'REPLACE')

      I_random_seed = int(2.d0**31 * grnd())
      if (I_random_seed.le.0) I_random_seed = (2**31 - 1) + I_random_seed

      WRITE (9, '(2x,i14,2x,i14)') T_cntr,           I_random_seed
      WRITE (9, '(2x,i14,2x,i14)') Start_diag_Tcntr, current_snap
      WRITE (9, '(2x,i14,2x,i14)') N_part(1),        N_part(2)
      WRITE (9, '(2x,i14,2x,i14)') Q_left,           Q_right
      WRITE (9, '(2x,e14.7,2x,e14.7,2x,e14.7)') full_Q_left,  full_Q_right , Q_ext

      WRITE (9, '(2x,i14,2x,i14)')     N_of_saved_records, text_output_counter
      WRITE (9, '(2x,e14.7,2x,e14.7)') Energy_full_eV,     Init_energy_full_eV
      WRITE (9, '(2x,e14.7,2x,e14.7)') Energy_pot_eV(1),   Energy_pot_eV(2)
      WRITE (9, '(2x,e14.7,2x,e14.7)') Energy_heat_eV(1),  Energy_heat_eV(2)
      WRITE (9, '(2x,e14.7,2x,e14.7)') Energy_kin_eV(1),   Energy_kin_eV(2)
      WRITE (9, '(2x,e14.7,2x,e14.7)') Energy_wall_eV(1),  Energy_wall_eV(2)
      WRITE (9, '(2x,e14.7,2x,e14.7)') Energy_emit_eV(1),  Energy_emit_eV(2)
      WRITE (9, '(2x,e14.7,2x,e14.7)') Energy_coll_eV(1),  Energy_coll_eV(2)
      WRITE (9, '(2x,e14.7,2x,e14.7)') prev_Q_left,      prev_Q_right
      CLOSE (9, STATUS = 'KEEP')

   ELSE

      OPEN (9, FILE = check_g_filename, STATUS = 'REPLACE')

      I_random_seed = int(2.d0**31 * grnd())
      if (I_random_seed.le.0) I_random_seed = (2**31 - 1) + I_random_seed

      WRITE (9, '(2x,i14,2x,i14)') T_cntr,           I_random_seed
      WRITE (9, '(2x,i14,2x,i14)') Start_diag_Tcntr, current_snap
      WRITE (9, '(2x,i14,2x,i14)') N_part(1),        N_part(2)
      WRITE (9, '(2x,i14,2x,i14)') Q_left,           Q_right
      CLOSE (9, STATUS = 'KEEP')

      OPEN (9, FILE = check_e_filename, STATUS = 'REPLACE')
      DO n = 1, N_part(1)
         WRITE (9, '(5(1x,e16.9),1x,i4)')  X_of_spec(1)%part(n), &
         & VX_of_spec(1)%part(n), &
         & VY_of_spec(1)%part(n), &
         & VZ_of_spec(1)%part(n), &
         & AX_of_spec(1)%part(n), &
         & Tag_of_spec(1)%part(n)
      END DO
      CLOSE (9, STATUS = 'KEEP')

      IF (N_spec.EQ.2) THEN
         OPEN (9, FILE = check_i_filename, STATUS = 'REPLACE')
         DO n = 1, N_part(2)
            WRITE (9, '(5(1x,e16.9),1x,i4)')  X_of_spec(2)%part(n), &
            & VX_of_spec(2)%part(n), &
            & VY_of_spec(2)%part(n), &
            & VZ_of_spec(2)%part(n), &
            & AX_of_spec(2)%part(n), &
            & Tag_of_spec(2)%part(n)
         END DO
         CLOSE (9, STATUS = 'KEEP')
      END IF

   END IF

   PRINT '(2x,"Process ",i3," : created the new checkpoint...")', Rank_of_process

   Save_check_Tcntr = Save_check_Tcntr + SaveCheck_step

END SUBROUTINE SAVE_CHECKPOINT
