!==================================================================================
subroutine FINAL_PUSH
!! Performs the final push of particles. This procedure follows equation 3.8 in Sydorenko's
!! thesis. The "streaming position" refers to the position of the particle after the
!! predicting push while the "temp position" refers to the position of the particle after
!! the final push, before the position is checked for collisions with the wall.
!!
!! The procedure calculates an acceleration that is used to update particles in the
!! PREDICTING_PUSH subroutine.

   use CurrentProblemValues
   use Diagnostics
   use SEEmission
   use ElectronInjection
   use BeamInPlasma

   implicit none

   integer s            ! species type
   integer k            ! particle

   real(8) x_before_prepush, left_node, right_node
   real(8) alfa_x, alfa_y, alfa_x2, invtheta
   real(8) KvEx, KvEy, KvEz, KxEx

! temporary values
   integer n_str_left     ! node to the left of the particle at the streaming position
   integer n_left         ! node to the left of the particle at the final position
   real(8) Xtmp           ! X-coordinate
   real(8) Xstr           ! streaming coordinate before correction is applied
   real(8) EXstr          ! streaming X-electric field
   integer nodes_crossed, j, node

   real(8) v_2            ! squared absolute value of velocity

! functions
   real(8) Bx_gauss
   real(8) By_gauss

   Q_strm_spec = 0.0_8  ! Clean the densities because we will need non-streaming densities
   ! in the near-wall nodes

! Cycle over electrons
   s = 1
   k = 0

   do while (k.LT.N_part(s))

      k = k + 1

      ! ### restore coordinate before the pre-push to get the magnetic field
      x_before_prepush = X_of_spec(s)%part(k) - KVx * VX_of_spec(s)%part(k)
      if (PeriodicBoundaryFlag.EQ.1) then
         if (x_before_prepush.GT.DBLE(N_cells)) then
            x_before_prepush = x_before_prepush - DBLE(N_cells)
         else if (x_before_prepush.LT.0.0_8) then
            x_before_prepush = x_before_prepush + DBLE(N_cells)
         end if
      end if

      ! ### prepare coefficients of dynamic equations depending on the magnetic field
      alfa_x = alfa_xy(s) * Bx_gauss(x_before_prepush)
      alfa_y = alfa_xy(s) * By_gauss(x_before_prepush)
      alfa_x2 = alfa_x**2
      invtheta = 1.0_8 / (1.0_8 + alfa_x2 + alfa_y**2)

      KvEx = KvE_xyz(s) * ((1.0_8 + alfa_x2) * invtheta) ! = KvE_xyz(s) if B_ext = 0
      KvEy = KvE_xyz(s) *   (alfa_x * alfa_y * invtheta) ! = 0.0_8 - " -
      KvEz = KvE_xyz(s) *            (alfa_y * invtheta) ! = 0.0_8 - " -

      KxEx = KvEx * factor_KxEx                          ! = KvE_xyz(s) * factor_KxEx if B_ext = 0
      ! ###

      if (X_of_spec(s)%part(k).GT.N_cells) then     !
         EXstr = EX(N_cells)                        ! for particles that crossed the plasma borders during the predicting push
      else if (X_of_spec(s)%part(k).LT.0.0_8) then  !
         EXstr = EX(0)                              !
      else
         n_str_left   = INT(X_of_spec(s)%part(k))           ! obtained after pre-push, corresponds to streaming position
         if (n_str_left.EQ.N_cells) n_str_left = N_cells - 1
         EXstr = EX(n_str_left) + GradEX(n_str_left) * (X_of_spec(s)%part(k) - n_str_left)
      end if

      ! velocity correction
      VX_of_spec(s)%part(k) = VX_of_spec(s)%part(k) + KvEx * EXstr
      VY_of_spec(s)%part(k) = VY_of_spec(s)%part(k) + KvEy * EXstr
      VZ_of_spec(s)%part(k) = VZ_of_spec(s)%part(k) + KvEz * EXstr

      if (UseSmartTagsFlag.EQ.1) then
         ! if tracing of whether a particle changed its direction of propagation along the x-axis is requested
         ! clear the tag of the particle if the x-velocity changes sign after x-acceleration (i.e. particle turns back or stops)
         if (prev_VX_of_spec(s)%part(k).GE.0.0_8) then                                          !
            if (VX_of_spec(s)%part(k).LE.0.0_8) Tag_of_spec(s)%part(k) = 0     !
         else                                                                  ! diagnostics
            if (VX_of_spec(s)%part(k).GE.0.0_8) Tag_of_spec(s)%part(k) = 0     !
         end if                                                                !
      end if

      ! x-move
      Xstr = X_of_spec(s)%part(k)
      Xtmp = X_of_spec(s)%part(k) + KxEx * EXstr ! *** Eqn 3.8 in Sydorenko's thesis ***
      X_of_spec(s)%part(k) = Xtmp ! before checks for wall collision, periodicity

      ! find number of nodes (with a sign) passed by the particle, not accounting for periodic b.c. yet:
      left_node = dble(int(Xtmp))
      right_node = 1.0_8 + left_node
      nodes_crossed = int(Xtmp) - int(x_before_prepush) !x_before_prepush defined to always be inside the domain
      if (nodes_crossed .ge. 1) then
         do j = 1, nodes_crossed, 1
            node = left_node - j + 1
            if (node.gt.0 .and. node.lt.N_cells) NVX_mesh(node, s) = NVX_mesh(node, s) + 1.0_8  ! add flux to each node crossed in positive direction
         end do
      end if

      if (nodes_crossed .le. -1) then
         do j = 1, -nodes_crossed, 1
            node = right_node + j - 1
            if (node.lt.N_cells .and. node.gt.0) NVX_mesh(node, s) = NVX_mesh(node, s) - 1.0_8 ! crossed in negative direction
         end do
      end if

      if (Xtmp.GT.N_cells) then
         ! if crossed the right boundary, electrons:
         v_2 = VX_of_spec(s)%part(k)**2 + VY_of_spec(s)%part(k)**2 + VZ_of_spec(s)%part(k)**2      !
         Rate_energy_rightwall(s) = Rate_energy_rightwall(s) + v_2                                 ! diagnostics
         Rate_number_rightwall(s) = Rate_number_rightwall(s) + 1                                   !
         NVX_mesh(N_cells, s) = 1. + NVX_mesh(N_cells, s)                                          !
         if (PeriodicBoundaryFlag.EQ.1) then
            X_of_spec(s)%part(k) = Xtmp - N_cells   ! move particle to the left end of the system
            Rate_energy_leftemit(s) = Rate_energy_leftemit(s) + v_2                                ! diagnostics
            Rate_number_leftemit(s) = Rate_number_leftemit(s) + 1                                  !

            if (Xtmp.GT.dble(1 + N_cells)) then                                                    !
               do j = 1 + N_cells, int(Xtmp)                                                       !
                  NVX_mesh(j - N_cells, s) = 1. + NVX_mesh(j - N_cells, s)                         !
               end do                                                                              !
            end if                                                                                 !
         else
            if (Xstr.LE.N_cells) then
               ! since we are here, collision occurred NOT at the pre-push, but after the final correction was applied
               ! so we must update the wall charge
               Q_right = Q_right + Qs(s)
            end if
            select case (Tag_of_spec(s)%part(k))                                                   !
             case (eTag_Emit_Left)                                                                 !
               prie_right_from_left_count   = prie_right_from_left_count   + 1                     !
               prie_right_from_left_energy  = prie_right_from_left_energy  + v_2                   ! diagnostics
             case (eTag_Coll_Neutral)                                                              !
               prie_right_after_coll_count  = prie_right_after_coll_count  + 1                     !
               prie_right_after_coll_energy = prie_right_after_coll_energy + v_2                   !
            end select                                                                             !
            call ADD_PRIMARY_E_TO_RIGHT_DF(VX_of_spec(s)%part(k), VY_of_spec(s)%part(k), VZ_of_spec(s)%part(k))      ! diagnostics
            call PROCESS_ELECTRON_COLLISION_WITH_WALL(s, Xtmp, k)
            call SUBSTITUTE_LEAVING_PARTICLE(s, k)
            k = k - 1 ! Move the index back one particle, to iterate the end particle now in the kth position
            N_part(s) = N_part(s) - 1
            cycle
         end if
      else if (Xtmp.LT.0.0_8) then
         ! if crossed the left boundary, electrons:
         v_2 = VX_of_spec(s)%part(k)**2 + VY_of_spec(s)%part(k)**2 + VZ_of_spec(s)%part(k)**2      !
         Rate_energy_leftwall(s)  = Rate_energy_leftwall(s)  + v_2                                 ! diagnostics
         Rate_number_leftwall(s)  = Rate_number_leftwall(s)  + 1                                   !
         NVX_mesh(0, s) = -1. + NVX_mesh(0, s)                                                     !
         if (PeriodicBoundaryFlag.EQ.1) then
            X_of_spec(s)%part(k) = N_cells + Xtmp    ! move particle to the right end of the system
            Rate_energy_rightemit(s) = Rate_energy_rightemit(s) + v_2                              ! diagnostics
            Rate_number_rightemit(s) = Rate_number_rightemit(s) + 1                                !

            if (Xtmp .lt. -1.0_8) then                                                             !
               do j = 1 + int(dble(N_cells) + Xtmp), N_cells - 1                                   !
                  NVX_mesh(j, s) = -1. + NVX_mesh(j, s)                                            !
               end do                                                                              !
            end if                                                                                 !
         else
            if (Xstr.GE.0.0_8) then
               ! since we are here, collision occurred NOT at the pre-push, but after the final correction was applied
               ! so we must update the wall charge
               Q_left = Q_left + Qs(s)
            end if
            select case (Tag_of_spec(s)%part(k))                                                   !
             case (eTag_Emit_Right)                                                                !
               prie_left_from_right_count  = prie_left_from_right_count  + 1                       !
               prie_left_from_right_energy = prie_left_from_right_energy + v_2                     ! diagnostics
             case (eTag_Coll_Neutral)                                                              !
               prie_left_after_coll_count  = prie_left_after_coll_count  + 1                       !
               prie_left_after_coll_energy = prie_left_after_coll_energy + v_2                     !
            end select                                                                             !
            call ADD_PRIMARY_E_TO_LEFT_DF(VX_of_spec(s)%part(k), VY_of_spec(s)%part(k), VZ_of_spec(s)%part(k))    ! diagnostics
            call PROCESS_ELECTRON_COLLISION_WITH_WALL(s, Xtmp, k)
            call SUBSTITUTE_LEAVING_PARTICLE(s, k)
            k = k - 1
            N_part(s) = N_part(s) - 1
            cycle
         end if
      else if (Xstr.GT.N_cells) then
         ! crossed the right boundary at the pre-push but after correction returned to the simulation area
         ! so we must remove contribution to the wall charge
         Q_right = Q_right - Qs(s)
      else if (Xstr.LT.0.0_8) then
         ! crossed the left boundary at the pre-push but after correction returned to the simulation area
         ! so we must remove contribution to the wall charge
         Q_left = Q_left - Qs(s)
      end if

      ! calculate the non-streaming density (not used for finding electric field?)
      n_left = INT(Xtmp)
      if (n_left.EQ.N_cells) n_left = n_left - 1    ! correct, if the particle is EXACTLY at the right wall
      Q_strm_spec(n_left,   s)  = Q_strm_spec(n_left, s) + (n_left + 1 - Xtmp)
      Q_strm_spec(n_left+1, s)  = Q_strm_spec(n_left+1, s) + (Xtmp - n_left)

      ! note that for periodic boundary conditions, accelerations of electrons crossing the boundaries are updated here as well
      AX_of_spec(s)%part(k) = 0.5_8 * (AX_of_spec(s)%part(k) + QMs(s) * EXstr) ! *** Eqn 3.8 in Sydorenko's thesis ***

   end do ! End of cycle over electrons

   ! Cycle over ions
   if (N_spec.ge.2) then

      s = 2
      k = 0
      do while (k.LT.N_part(s))

         k = k + 1

         ! ### restore coordinate before the pre-push to get the magnetic field, ions:
         x_before_prepush = X_of_spec(s)%part(k) - KVx * VX_of_spec(s)%part(k)
         if (PeriodicBoundaryFlag.EQ.1) then
            if (x_before_prepush.GT.DBLE(N_cells)) then
               x_before_prepush = x_before_prepush - DBLE(N_cells)
            else if (x_before_prepush.LT.0.0_8) then
               x_before_prepush = x_before_prepush + DBLE(N_cells)
            end if
         end if

         ! ### prepare coefficients of dynamic equations depending on the magnetic field
         alfa_x = alfa_xy(s) * Bx_gauss(x_before_prepush)
         alfa_y = alfa_xy(s) * By_gauss(x_before_prepush)
         alfa_x2 = alfa_x**2
         invtheta = 1.0_8 / (1.0_8 + alfa_x2 + alfa_y**2)

         KvEx = KvE_xyz(s) * ((1.0_8 + alfa_x2) * invtheta)
         KvEy = KvE_xyz(s) *   (alfa_x * alfa_y * invtheta)
         KvEz = KvE_xyz(s) *            (alfa_y * invtheta)

         KxEx = KvEx * factor_KxEx
         ! ###

         if (X_of_spec(s)%part(k).GT.N_cells) then              !
            EXstr = EX(N_cells)                                 ! for particles that crossed the
         else if (X_of_spec(s)%part(k).LT.0.0_8) then           ! plasma borders during the predicting push
            EXstr = EX(0)                                       !
         else
            n_str_left   = INT(X_of_spec(s)%part(k))            ! obtained after pre-push, corresponds to streaming position
            if (n_str_left.EQ.N_cells) n_str_left = N_cells - 1
            EXstr        = EX(n_str_left) + GradEX(n_str_left) * (X_of_spec(s)%part(k) - n_str_left)
         end if

         ! velocity correction
         VX_of_spec(s)%part(k) = VX_of_spec(s)%part(k) + KvEx * EXstr
         VY_of_spec(s)%part(k) = VY_of_spec(s)%part(k) + KvEy * EXstr
         VZ_of_spec(s)%part(k) = VZ_of_spec(s)%part(k) + KvEz * EXstr

         ! x-move
         Xstr = X_of_spec(s)%part(k)
         Xtmp = X_of_spec(s)%part(k) + KxEx * EXstr
         X_of_spec(s)%part(k) = Xtmp

         ! find number of nodes (with a sign) passed by the particle, not accounting for periodic b.c. yet:
         left_node = dble(int(Xtmp))
         right_node = 1.0_8 + left_node
         nodes_crossed = int(Xtmp) - int(x_before_prepush)

         if (nodes_crossed.ge. 1) then
            do j = 1, nodes_crossed
               node = left_node - j + 1
               if (node.gt.0 .and. node.lt.N_cells) NVX_mesh(node,s) = NVX_mesh(node,s) + 1.  ! crossed in positive direction
            end do
         end if

         if (nodes_crossed .le. -1) then
            do j = 1, -nodes_crossed
               node = right_node + j - 1
               if (node.lt.N_cells .and. node.gt.0) NVX_mesh(node,s) = NVX_mesh(node,s) - 1. ! crossed in negative direction
            end do
         end if

         if (Xtmp.GT.N_cells) then
            ! if crossed the right boundary, ions:
            v_2 = VX_of_spec(s)%part(k)**2 + VY_of_spec(s)%part(k)**2 + VZ_of_spec(s)%part(k)**2   !
            Rate_energy_rightwall(s) = Rate_energy_rightwall(s) + v_2                              ! diagnostics
            Rate_number_rightwall(s) = Rate_number_rightwall(s) + 1                                !
            NVX_mesh(N_cells, s) = 1. + NVX_mesh(N_cells, s)
            if (PeriodicBoundaryFlag.EQ.1) then
               X_of_spec(s)%part(k) = Xtmp - N_cells   ! move particle to the left end of the system
               Rate_energy_leftemit(s) = Rate_energy_leftemit(s) + v_2                             ! diagnostics
               Rate_number_leftemit(s) = Rate_number_leftemit(s) + 1                               !

               if (Xtmp .gt. dble(1 + N_cells)) then                                               !
                  do j = 1 + N_cells, int(Xtmp)                                                    !
                     NVX_mesh(j - N_cells, s) = 1. + NVX_mesh(j - N_cells, s)                      !
                  end do                                                                           !
               else                                                                                !
               end if                                                                              !
            else
               if (Xstr.LE.N_cells) then
                  ! since we are here, collision occurred NOT at the pre-push, but after the final correction was applied
                  ! so we must update the wall charge
                  Q_right = Q_right + Qs(s)
               end if
               call PROCESS_ION_COLLISION_WITH_RIGHT_WALL(Xtmp, VX_of_spec(s)%part(k), &
               & VY_of_spec(s)%part(k), VZ_of_spec(s)%part(k), v_2, s)
               call SUBSTITUTE_LEAVING_PARTICLE(s, k)
               k = k - 1
               N_part(s) = N_part(s) - 1
               cycle
            end if
         else if (Xtmp.LT.0.0_8) then
            ! if crossed the left boundary, ions:
            v_2 = VX_of_spec(s)%part(k)**2 + VY_of_spec(s)%part(k)**2 + VZ_of_spec(s)%part(k)**2   !
            Rate_energy_leftwall(s) = Rate_energy_leftwall(s) + v_2                                ! diagnostics
            Rate_number_leftwall(s) = Rate_number_leftwall(s) + 1                                  !
            NVX_mesh(0, s) = -1. + NVX_mesh(0, s)
            if (PeriodicBoundaryFlag.EQ.1) then
               X_of_spec(s)%part(k) = N_cells + Xtmp    ! move particle to the right end of the system
               Rate_energy_rightemit(s) = Rate_energy_rightemit(s) + v_2                           ! diagnostics
               Rate_number_rightemit(s) = Rate_number_rightemit(s) + 1                             !
               if (Xtmp .lt. -1.0_8) then                                                          !
                  do j = 1 + int( dble(N_cells) + Xtmp ), N_cells - 1                              !
                     NVX_mesh(j, s) = -1. + NVX_mesh(j, s)                                         !
                  end do                                                                           !
               else                                                                                !
               end if                                                                              !
            else
               if (Xstr.GE.0.0_8) then
                  ! since we are here, collision occurred NOT at the pre-push, but after the final correction was applied
                  ! so we must update the wall charge
                  Q_left = Q_left + Qs(s)  ! only if it was not updated at the pre-push
               end if
               call PROCESS_ION_COLLISION_WITH_LEFT_WALL(Xtmp, VX_of_spec(s)%part(k), &
               & VY_of_spec(s)%part(k), VZ_of_spec(s)%part(k), v_2, s)
               call SUBSTITUTE_LEAVING_PARTICLE(s, k)
               k = k - 1
               N_part(s) = N_part(s) - 1
               cycle
            end if
         else if (Xstr.GT.N_cells) then
            ! crossed the right boundary at the pre-push but after correction returned to the simulation area
            ! so we must remove contribution to the wall charge
            Q_right = Q_right - Qs(s)
         else if (Xstr.LT.0.0_8) then
            ! crossed the left boundary at the pre-push but after correction returned to the simulation area
            ! so we must remove contribution to the wall charge
            Q_left = Q_left - Qs(s)
         end if

         ! calculate the non-streaming density
         n_left = INT(Xtmp)
         if (n_left.EQ.N_cells) n_left = n_left - 1    ! correct, if the particle is EXACTLY at the right wall
         Q_strm_spec(n_left,   s)  = Q_strm_spec(n_left, s) + (n_left + 1 - Xtmp)
         Q_strm_spec(n_left+1, s)  = Q_strm_spec(n_left+1, s) + (Xtmp - n_left)

         ! note that for periodic boundary conditions, accelerations of ions crossing the boundarise are updated here as well
         AX_of_spec(s)%part(k) = 0.5_8 * (AX_of_spec(s)%part(k) + QMs(s) * EXstr) ! *** Eqn 3.8 in Sydorenko's thesis ***

      end do ! End of cycle over ions

   end if

   ! this should not be necessary, just a paranoidal precaution...
   if (PeriodicBoundaryFlag.EQ.1) then
      ! for periodic boundaries prevent charge accumulation
      Q_right = 0
      Q_left = 0
   end if

end subroutine FINAL_PUSH

!===============================================================================
subroutine PREDICTING_PUSH
!! Performs the predicting push of particles. This procedure follows equation 3.7 in
!! Sydorenko's thesis. The "streaming position" refers to the position of the particle
!! after the predicting push. Wall charge is updated here for Poisson's equation calculation,
!! but particles are not removed or undergo SEE until the final push.
!!
!! The procedure uses an acceleration that is calculated in the FINAL_PUSH subroutine.
!! As such, this should be called after FINAL_PUSH.

   use CurrentProblemValues
   use Diagnostics
   use SEEmission
   use ElectronInjection
   use BeamInPlasma

   implicit none

   integer s            ! species type
   integer k            ! particle

   real(8) alfa_x, alfa_y, alfa_x2, alfa_y2, theta2, invtheta
   real(8) K11, K12, K13, K22, K23, K33 ! Coefficients (before V^{n-1/2}) for equations of motion
                                        ! (see Eqn 3.2 in Sydorenko's thesis)
   real(8) A11, A21, A31 ! Coefficients before Ax^{n-1} for equations of motion
   real(8) A13, A23, A33 ! Coefficients with Ez_external for equations of motion

   real(8) new_Vx, new_Vy

   integer n_left, n_right            ! bounding nodes

   real(8) Xstr, delta

   ! functions
   real(8) Bx_gauss
   real(8) By_gauss

   ! if necessary set flags for accumulation of distribution functions at the walls
   call REQUEST_WALL_DF_FOR_SNAPSHOT

   N_inject = 0         ! Clean the counters
   Q_strm_spec = 0.0_8  ! Clean the densities

   ! clean wall charge density accumulators if we simulate a situation with given wall potential or external circuit
   ! because in this case the final wall charge density at the end of previous time step is calculated from the Gauss law
   if (BC_flag.eq.0 .or. BC_flag.eq.2) then
      Q_left = 0
      Q_right = 0
   end if

   ! periodic boundary ---------------------------------------
   if (PeriodicBoundaryFlag.EQ.1) then

      ! cycle over species
      do s = 1, N_spec

         ! cycle over particles
         do k = 1, N_part(s)

            ! ### prepare coefficients of dynamic equations depending on the magnetic field
            alfa_x = alfa_xy(s) * Bx_gauss(X_of_spec(s)%part(k))
            alfa_y = alfa_xy(s) * By_gauss(X_of_spec(s)%part(k))
            alfa_x2 = alfa_x**2
            alfa_y2 = alfa_y**2
            theta2 = alfa_x2 + alfa_y2
            invtheta = 1.0_8 / (1.0_8 + theta2)

            K11 = (1.0_8 - theta2 + 2.0_8 * alfa_x2) * invtheta ! = 1.0_8 if B_ext = 0
            K12 =            2.0_8 * alfa_x * alfa_y * invtheta ! = 0.0_8 - " -
            K13 =                     2.0_8 * alfa_y * invtheta ! = 0.0_8 - " -
            K22 = (1.0_8 - theta2 + 2.0_8 * alfa_y2) * invtheta ! = 1.0_8 - " -
            K23 =                     2.0_8 * alfa_x * invtheta ! = 0.0_8 - " -
            K33 =                   (1.0_8 - theta2) * invtheta ! = 1.0_8 - " -

            A11 = A_123_1 * ((1.0_8 + alfa_x2) * invtheta)      ! = A_123_1 if B_ext = 0
            A21 = A_123_1 *   (alfa_x * alfa_y * invtheta)      ! = 0.0_8   - " -
            A31 = A_123_1 *            (alfa_y * invtheta)      ! = 0.0_8   - " -

            A13 = A_123_3(s) * (alfa_y * invtheta)              ! = 0.0_8 if E_ext = 0 or B_ext = 0
            A23 = A_123_3(s) * (alfa_x * invtheta)              ! = 0.0_8 - " -
            A33 = A_123_3(s)           * invtheta               ! = A_123_3(s) if B_ext = 0 and = 0.0_8 if E_ext = 0
            ! ###

            ! *** Pushes follow equation 3.7 in Sydorenko's thesis ***
            ! x-acceleration
            new_Vx = K11 * VX_of_spec(s)%part(k) + &            ! = VX_of_spec(s)%part(k)             if B_ext = 0
            & K12 * VY_of_spec(s)%part(k) - &                   ! = + 0.0_8                           - " -
            & K13 * VZ_of_spec(s)%part(k) + &                   ! = - 0.0_8                           - " -
            & A11 * AX_of_spec(s)%part(k) - &                   ! = + A_123_1 * AX_of_spec(s)%part(k) - " -
            & A13                                               ! = - 0.0_8                           if E_ext = 0

            ! y-acceleration
            new_Vy = K12 * VX_of_spec(s)%part(k) + &
            & K22 * VY_of_spec(s)%part(k) + &
            & K23 * VZ_of_spec(s)%part(k) + &
            & A21 * AX_of_spec(s)%part(k) + &
            & A23

            ! z-acceleration
            VZ_of_spec(s)%part(k) = K13 * VX_of_spec(s)%part(k) - &
            & K23 * VY_of_spec(s)%part(k) + &
            & K33 * VZ_of_spec(s)%part(k) + &
            & A31 * AX_of_spec(s)%part(k) + &
            & A33

            prev_VX_of_spec(s)%part(k) = VX_of_spec(s)%part(k)
            VX_of_spec(s)%part(k) = new_Vx
            VY_of_spec(s)%part(k) = new_Vy

            ! x-move
            Xstr = X_of_spec(s)%part(k) + KVx * VX_of_spec(s)%part(k) ! *** Eqn 3.7 in Sydorenko's thesis ***

            if (Xstr.GT.N_cells) then
               ! crossed the right boundary
               Xstr = Xstr - N_cells   ! move particle to the left end of the system
            else if (Xstr.LT.0.0_8) then
               ! crossed the left boundary
               Xstr = Xstr + N_cells   ! move particle to the right end of the system
            end if
            X_of_spec(s)%part(k) = Xstr

            n_left = INT(Xstr)
            if (n_left.EQ.N_cells) n_left = n_left - 1    ! correct, if the particle is EXACTLY at the right wall
            n_right = n_left + 1

            Q_strm_spec(n_left, s)  = Q_strm_spec(n_left, s)  + n_right - Xstr
            Q_strm_spec(n_right, s) = Q_strm_spec(n_right, s) + Xstr - n_left
         end do
      end do

      ! non-periodic boundary ------------------------------------------
   else
      ! modified 10/30/2020

      ! cycle over species
      do s = 1, N_spec

         ! cycle over particles
         do k = 1, N_part(s)

            ! ### prepare coefficients of dynamic equations depending on the magnetic field
            alfa_x = alfa_xy(s) * Bx_gauss(X_of_spec(s)%part(k))
            alfa_y = alfa_xy(s) * By_gauss(X_of_spec(s)%part(k))
            alfa_x2 = alfa_x**2
            alfa_y2 = alfa_y**2
            theta2 = alfa_x2 + alfa_y2
            invtheta = 1.0_8 / (1.0_8 + theta2)

            K11 = (1.0_8 - theta2 + 2.0_8 * alfa_x2) * invtheta ! = 1.0_8 if B_ext = 0
            K12 =            2.0_8 * alfa_x * alfa_y * invtheta ! = 0.0_8 - " -
            K13 =                     2.0_8 * alfa_y * invtheta ! = 0.0_8 - " -
            K22 = (1.0_8 - theta2 + 2.0_8 * alfa_y2) * invtheta ! = 1.0_8 - " -
            K23 =                     2.0_8 * alfa_x * invtheta ! = 0.0_8 - " -
            K33 =                   (1.0_8 - theta2) * invtheta ! = 1.0_8 - " -

            A11 = A_123_1 * ((1.0_8 + alfa_x2) * invtheta)      ! = A_123_1 if B_ext = 0
            A21 = A_123_1 *   (alfa_x * alfa_y * invtheta)      ! = 0.0_8   - " -
            A31 = A_123_1 *            (alfa_y * invtheta)      ! = 0.0_8   - " -

            A13 = A_123_3(s) * (alfa_y * invtheta)              ! = 0.0_8 if E_ext = 0 or B_ext = 0
            A23 = A_123_3(s) * (alfa_x * invtheta)              ! = 0.0_8 - " -
            A33 = A_123_3(s)           * invtheta               ! = A_123_3(s) if B_ext = 0 and = 0.0_8 if E_ext = 0
            ! ###

            ! *** Pushes follow equation 3.7 in Sydorenko's thesis ***
            ! x-acceleration
            new_Vx = K11 * VX_of_spec(s)%part(k) + &            ! = VX_of_spec(s)%part(k)             if B_ext = 0
            & K12 * VY_of_spec(s)%part(k) - &                   ! = + 0.0_8                           - " -
            & K13 * VZ_of_spec(s)%part(k) + &                   ! = - 0.0_8                           - " -
            & A11 * AX_of_spec(s)%part(k) - &                   ! = + A_123_1 * AX_of_spec(s)%part(k) - " -
            & A13                                               ! = - 0.0_8                           if E_ext = 0

            ! y-acceleration
            new_Vy = K12 * VX_of_spec(s)%part(k) + &
            & K22 * VY_of_spec(s)%part(k) + &
            & K23 * VZ_of_spec(s)%part(k) + &
            & A21 * AX_of_spec(s)%part(k) + &
            & A23

            ! z-acceleration
            VZ_of_spec(s)%part(k) = K13 * VX_of_spec(s)%part(k) - &
            & K23 * VY_of_spec(s)%part(k) + &
            & K33 * VZ_of_spec(s)%part(k) + &
            & A31 * AX_of_spec(s)%part(k) + &
            & A33

            prev_VX_of_spec(s)%part(k) = VX_of_spec(s)%part(k)
            VX_of_spec(s)%part(k) = new_Vx
            VY_of_spec(s)%part(k) = new_Vy

            ! x-move
            Xstr = X_of_spec(s)%part(k) + KVx * VX_of_spec(s)%part(k) ! *** Eqn 3.7 in Sydorenko's thesis ***
            X_of_spec(s)%part(k) = Xstr

            if (Xstr.GT.N_cells) then
               ! crossed the right boundary
               Q_right = Q_right + Qs(s)
            else if (Xstr.LT.0.0_8) then
               ! crossed the left boundary
               Q_left = Q_left + Qs(s)
            else
               ! stays inside the system
               n_left = INT(Xstr)
               if (n_left.EQ.N_cells) n_left = n_left - 1    ! correct, if the particle is EXACTLY at the right wall
               n_right = n_left + 1
               ! modification 10/30/2020:
               delta = Xstr - n_left - 0.5_8
               if (delta.ge.0.) then
                  Q_strm_spec(n_left, s)  = Q_strm_spec(n_left, s) + 0.5 * (1.- delta)**2
                  Q_strm_spec(n_right, s) = Q_strm_spec(n_right, s) + 1. - 0.5 * (1. - delta)**2 - 0.5 * delta**2
                  if (n_left.le.N_cells - 2) then
                     Q_strm_spec(n_right + 1, s) = Q_strm_spec(n_right + 1, s) + 0.5 * delta ** 2
                  end if
               else
                  Q_strm_spec(n_right, s) = Q_strm_spec(n_right, s) + 0.5 * (1. + delta)**2
                  Q_strm_spec(n_left, s)  = Q_strm_spec(n_left, s) + 1. - 0.5 * (1. + delta)**2 - 0.5 * delta**2
                  if (n_left.ge.1) then
                     Q_strm_spec(n_left - 1, s)  = Q_strm_spec(n_left - 1, s) + 0.5 * delta**2
                  end if
               end if
               !changed from:
               !!              Q_strm_spec(n_left, s)  = Q_strm_spec(n_left, s)  + n_right - Xstr
               !!              Q_strm_spec(n_right, s) = Q_strm_spec(n_right, s) + Xstr - n_left
            end if
         end do
      end do

   end if

end subroutine PREDICTING_PUSH

!===============================================================================
subroutine SUBSTITUTE_LEAVING_PARTICLE(s, k)
!! Removes the species s particle with index k. If this is not the last particle in
!! the species list, the index k particle is replaced by the last particle in the array
!! and the last particle is removed.

   use CurrentProblemValues, only : N_part, X_of_spec, VX_of_spec, VY_of_spec, VZ_of_spec, AX_of_spec, Tag_of_spec, prev_VX_of_spec
   implicit none

   integer s            ! species type
   integer k            ! index of particle to be removed

   if (k.LT.N_part(s)) then ! if the data are not the last active one
      X_of_spec(s)%part(k)       = X_of_spec(s)%part(N_part(s))        !
      VX_of_spec(s)%part(k)      = VX_of_spec(s)%part(N_part(s))       ! SET data equal
      VY_of_spec(s)%part(k)      = VY_of_spec(s)%part(N_part(s))       ! to data from LAST
      VZ_of_spec(s)%part(k)      = VZ_of_spec(s)%part(N_part(s))       ! ACTIVE position
      AX_of_spec(s)%part(k)      = AX_of_spec(s)%part(N_part(s))       !
      Tag_of_spec(s)%part(k)     = Tag_of_spec(s)%part(N_part(s))      !
      prev_VX_of_spec(s)%part(k) = prev_VX_of_spec(s)%part(N_part(s))  !
   end if                                                !

   X_of_spec(s)%part(N_part(s))       = 0.0_8 !
   VX_of_spec(s)%part(N_part(s))      = 0.0_8 ! CLEAR LAST
   VY_of_spec(s)%part(N_part(s))      = 0.0_8 ! ACTIVE data
   VZ_of_spec(s)%part(N_part(s))      = 0.0_8 !
   AX_of_spec(s)%part(N_part(s))      = 0.0_8 !
   Tag_of_spec(s)%part(N_part(s))     = 0     !
   prev_VX_of_spec(s)%part(N_part(s)) = 0.0_8 !

end subroutine SUBSTITUTE_LEAVING_PARTICLE

!===============================================================================
subroutine MANAGE_WALL_CHARGE_DENSITIES

   use mpi
   use ParallelOperationValues
   use CurrentProblemValues
   implicit none

!  INCLUDE 'mpif.h'

   real(8) rbufer(1:2), rbufer2(1:2)
   integer ierr

   if (BC_flag.ne.0 .and. BC_flag.ne.2 ) return

   rbufer  = 0.0_8
   rbufer2 = 0.0_8

   ! send volume charge densities at the walls to server
   if (Rank_of_process.GT.0) then

      rbufer(1) = Qs(1) * Q_strm_spec(0, 1)       + Qs(2) * Q_strm_spec(0, 2)
      rbufer(2) = Qs(1) * Q_strm_spec(N_cells, 1) + Qs(2) * Q_strm_spec(N_cells, 2)

      call MPI_REDUCE(rbufer, rbufer2, 2, MPI_doUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

   ! receive from clients volume charge densities assigned to wall nodes
   else

      call MPI_REDUCE(rbufer2, rbufer, 2, MPI_doUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      full_Q_left  = (F(0) - F(1)) / K_Q - rbufer(1)
      full_Q_right = (F(N_cells) - F(N_cells-1)) / K_Q - rbufer(2)

   end if

end subroutine MANAGE_WALL_CHARGE_DENSITIES
