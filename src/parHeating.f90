!----------------------------------
module heating_diagnostics

   real(8), allocatable :: joule_heat_Wm3(:,:)     !! Rate of Joule heating for each species across the domain at the most
                                                   !! recent timestep. Allocated to size (1:N_cells, 1:N_spec)
   real(8), allocatable :: joule_heat_Jm3(:,:)     !! Joule heat into each species across the domain at the prior time step.
                                                   !! Allocated to size (1:N_cells, 1:N_spec)
   integer :: heating_diag_start                   !! Start time for heat measurement
end module heating_diagnostics

!===================================================================================================
subroutine calc_joule_heating()
   !! Calculates the joule heat into electrons for each species on all clients
   use heating_diagnostics
   use CurrentProblemValues, only : N_cells, N_spec, N_part, T_cntr, X_of_spec, VX_of_spec, &
      EX, GradEX, E_scl_Vm, V_scl_ms, e_Cl, Qs!, VZ_of_spec, E_z_ext_Vm
   use ParallelOperationValues
   implicit none

   integer :: s, k, left_node
   real(8) :: x, vx!, vz
   real(8) :: heating   !! heating to a single particle
   real(8) :: EX_k      !! interpolated electric field at the particle

   ! Check if we should be here
   if (T_cntr.lt.heating_diag_start) return     ! Do not calculate heating before the start time
   if (Rank_of_process.eq.0) return             ! Accumulate data only if the process is a client

   ! Clear the heating array from the prior timestep
   joule_heat_Wm3 = 0.0_8

   ! Loop over species
   do s = 1, N_spec
      ! Loop over particles
      do k = 1, N_part(s)
         ! Initiate the temporary variables
         x            = X_of_spec(s)%part(k)
         left_node    = int(x)
         if (left_node.EQ.N_cells) then         ! If particle is at the right node exactly (but not collided yet)
            left_node = left_node - 1
         end if
         vx           = VX_of_spec(s)%part(k)
         ! vz           = VZ_of_spec(s)%part(k) ! Uncomment to account for an external E_z field

         ! Calculate the electric field at the particle
         if (x.GT.N_cells) then                        !
            EX_k = EX(N_cells)                         ! for particles that crossed the plasma borders during the predicting push
         else if (x.LT.0.0_8) then                     !
            EX_k = EX(0)                               !
         else
            EX_k = EX(left_node) + GradEX(left_node) * (x - left_node)
         end if

         ! Calculate dimensionless joule heating
         heating = EX_k * vx !+ vz * E_z_ext_Vm
         joule_heat_Wm3(left_node,s) = joule_heat_Wm3(left_node,s) + heating
      end do
      ! Dimensionalize
      joule_heat_Wm3(:,s) = joule_heat_Wm3(:,s) * E_scl_Vm * V_scl_ms * Qs(s) * e_Cl
   end do

end subroutine calc_joule_heating

!===================================================================================================
subroutine allocate_heating()
   !! Allocate the allocated arrays in the heating_diagnostics module
   use heating_diagnostics
   use CurrentProblemValues, only : N_cells, N_spec
   implicit none
   integer :: ALLOC_ERR

   ! Allocate the heating array
   if (.not.allocated(joule_heat_Wm3)) then
      allocate(joule_heat_Wm3(1:N_cells,1:N_spec), stat=ALLOC_ERR)
      if(ALLOC_ERR.ne.0) then
         print *, 'Error in ALLOCATE joule_heat_Wm3 !!!'
         print *, 'The program will be terminated now :('
         stop
      end if
   end if

   ! Allocate the heat array
   if (.not.allocated(joule_heat_Jm3)) then
      allocate(joule_heat_Jm3(1:N_cells,1:N_spec), stat=ALLOC_ERR)
      if(ALLOC_ERR.ne.0) then
         print *, 'Error in ALLOCATE joule_heat_Jm3 !!!'
         print *, 'The program will be terminated now :('
         stop
      end if
   end if

end subroutine allocate_heating

!===================================================================================================
subroutine deallocate_heating()
   !! Deallocate the allocated arrays in the heating_diagnostics module
   use heating_diagnostics
   implicit none
   integer :: DEALLOC_ERR

   ! Deallocate the heating array
   if (allocated(joule_heat_Wm3)) then
      deallocate(joule_heat_Wm3, stat=DEALLOC_ERR)
      if(DEALLOC_ERR.ne.0) then
         print *, 'Error in DEALLOCATE joule_heat_Wm3 !!!'
         print *, 'The program will be terminated now :('
         stop
      end if
   end if

   ! Deallocate the heating array
   if (allocated(joule_heat_Jm3)) then
      deallocate(joule_heat_Jm3, stat=DEALLOC_ERR)
      if(DEALLOC_ERR.ne.0) then
         print *, 'Error in DEALLOCATE joule_heat_Jm3 !!!'
         print *, 'The program will be terminated now :('
         stop
      end if
   end if

end subroutine deallocate_heating
