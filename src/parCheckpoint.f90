subroutine SAVE_CHECKPOINT

   use CurrentProblemValues
   use ParallelOperationValues
   use Diagnostics
   use Snapshots
   use mt19937

   implicit none

   integer n ! number of particle

   real R

   if (T_cntr.NE.Save_check_Tcntr) return

   if (Rank_of_process.EQ.0) then

      open (9, file = check_g_filename, status = 'replace')

      I_random_seed = int(2.d0**31 * grnd())
      if (I_random_seed.le.0) I_random_seed = 2147483647 + I_random_seed ! changed from = (2**31 - 1) + I_random_seed

      write (9, '(2x,i14,2x,i14)') T_cntr,           I_random_seed
      write (9, '(2x,i14,2x,i14)') Start_diag_Tcntr, current_snap
      write (9, '(2x,i14,2x,i14)') N_part(1),        N_part(2)
      write (9, '(2x,i14,2x,i14)') Q_left,           Q_right
      write (9, '(2x,e14.7,2x,e14.7,2x,e14.7)') full_Q_left,  full_Q_right , Q_ext

      write (9, '(2x,i14,2x,i14)')     N_of_saved_records, text_output_counter
      write (9, '(2x,e14.7,2x,e14.7)') Energy_full_eV,     Init_energy_full_eV
      write (9, '(2x,e14.7,2x,e14.7)') Energy_pot_eV(1),   Energy_pot_eV(2)
      write (9, '(2x,e14.7,2x,e14.7)') Energy_heat_eV(1),  Energy_heat_eV(2)
      write (9, '(2x,e14.7,2x,e14.7)') Energy_kin_eV(1),   Energy_kin_eV(2)
      write (9, '(2x,e14.7,2x,e14.7)') Energy_wall_eV(1),  Energy_wall_eV(2)
      write (9, '(2x,e14.7,2x,e14.7)') Energy_emit_eV(1),  Energy_emit_eV(2)
      write (9, '(2x,e14.7,2x,e14.7)') Energy_coll_eV(1),  Energy_coll_eV(2)
      write (9, '(2x,e14.7,2x,e14.7)') prev_Q_left,      prev_Q_right
      close (9, status = 'keep')

   ELSE

      open (9, file = check_g_filename, status = 'replace')

      I_random_seed = int(2.d0**31 * grnd())
      if (I_random_seed.le.0) I_random_seed = 2147483647 + I_random_seed ! changed from = (2**31 - 1) + I_random_seed

      write (9, '(2x,i14,2x,i14)') T_cntr,           I_random_seed
      write (9, '(2x,i14,2x,i14)') Start_diag_Tcntr, current_snap
      write (9, '(2x,i14,2x,i14)') N_part(1),        N_part(2)
      write (9, '(2x,i14,2x,i14)') Q_left,           Q_right
      close (9, status = 'keep')

      open (9, file = check_e_filename, status = 'replace')
      do n = 1, N_part(1)
         write (9, '(5(1x,e16.9),1x,i4)')  X_of_spec(1)%part(n), &
         & VX_of_spec(1)%part(n), &
         & VY_of_spec(1)%part(n), &
         & VZ_of_spec(1)%part(n), &
         & AX_of_spec(1)%part(n), &
         & Tag_of_spec(1)%part(n)
      end do
      close (9, status = 'keep')

      if (N_spec.EQ.2) then
         open (9, file = check_i_filename, status = 'replace')
         do n = 1, N_part(2)
            write (9, '(5(1x,e16.9),1x,i4)')  X_of_spec(2)%part(n), &
            & VX_of_spec(2)%part(n), &
            & VY_of_spec(2)%part(n), &
            & VZ_of_spec(2)%part(n), &
            & AX_of_spec(2)%part(n), &
            & Tag_of_spec(2)%part(n)
         end do
         close (9, status = 'keep')
      end if

   end if

   print '(2x,"Process ",i3," : created the new checkpoint...")', Rank_of_process

   Save_check_Tcntr = Save_check_Tcntr + SaveCheck_step

end subroutine SAVE_CHECKPOINT
