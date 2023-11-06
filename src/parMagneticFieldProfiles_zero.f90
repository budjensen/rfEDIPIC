!-------------------------------------------------
!
! UNIFORM BX
!
! BY IS PROPORTIONAL TO 1/x or exp(x/L) or to x^3
!
!-------------------------------------------------
!
MODULE MagneticFieldProfiles

! parameters for BX-profile
   REAL(8), PARAMETER :: Bx_gauss_0 = 0.0_8

! parameters for BY-profile
   REAL(8), PARAMETER :: By_gauss_0 = 10.0_8     ! wall x=0 !currently not used, By(x) defined as a cubic finction.
   REAL(8), PARAMETER :: By_gauss_N = 0.0_8    ! wall x=N_cells
!  integer, parameter :: N_0 = 50
   real(8), parameter :: x_0_m = 0.02 !no ramp if x_0_m exceeds discharge gap
   REAL(8) delta, b  ! By=b/(x+delta)
   REAL(8) alpha, x_0, x_max     ! By=By_gauss_N*exp(-alpha*(N_cells-x))


END MODULE MagneticFieldProfiles

!-------------------------------------------------
! this subroutine MUST be called before the first use of functions Bx_gauss and By_gauss
!
SUBROUTINE PREPARE_CALCULATION_OF_MAG_FIELD_PROFILES

   USE MagneticFieldProfiles
   USE CurrentProblemValues, ONLY : N_cells, T_cntr

   IMPLICIT NONE

! ### prepare coefficients for Bx_gauss

!    nothing to do here because Bx is uniform, change if necessary

! ### prepare coefficients for By_gauss

!  delta = DBLE(N_cells) * By_gauss_N / (By_gauss_0 - By_gauss_N)

!  b = delta * By_gauss_0

!  alpha = dlog(By_gauss_N/By_gauss_0)/dble(N_cells)

END SUBROUTINE PREPARE_CALCULATION_OF_MAG_FIELD_PROFILES


!-------------------------------------------------
! the result must be Bx-magnetic field in Gauss
!
! MODIFY THIS FUNCTION ACCORDING TO THE REQUIRED PROFILE OF BX
! AND, IF NECESSARY,
! MAKE APPROPRIATE CHANGES IN MODULE MagneticFieldProfiles AND
! SUBROUTINE PREPARE_CALCULATION_OF_MAG_FIELD_PROFILES
!
REAL(8) FUNCTION Bx_gauss(x)

   USE CurrentProblemValues, ONLY : N_cells, T_cntr
   USE MagneticFieldProfiles

   IMPLICIT NONE

   REAL(8) x   ! coordinate, dimensionless (in units of cell size)

! below a uniform Bx is implemented

   Bx_gauss = Bx_gauss_0

END FUNCTION Bx_gauss

!-------------------------------------------------
! the result must be By-magnetic field in Gauss
!
! MODIFY THIS FUNCTION ACCORDING TO THE REQUIRED PROFILE OF BY
! AND, IF NECESSARY,
! MAKE APPROPRIATE CHANGES IN MODULE MagneticFieldProfiles AND
! SUBROUTINE PREPARE_CALCULATION_OF_MAG_FIELD_PROFILES
!
REAL(8) FUNCTION By_gauss(x)

   USE CurrentProblemValues, ONLY : N_cells, T_cntr, delta_x_m
   USE MagneticFieldProfiles

   IMPLICIT NONE

   REAL(8) x   ! coordinate, dimensionless (in units of cell size)
!  By_gauss = By_gauss_N * exp(-alpha*(N_cells-x))
!  By_gauss = b / (x + delta)
!  By_gauss = By_gauss_0
!  x_0 = dble(N_cells - N_0)

   x_0 = x_0_m / delta_x_m
   x_max = dble(N_cells)
   if (x_0.gt.x_max) x_0 = x_max

   if (x.ge.x_0 .and. x_max.gt.x_0) then
      By_gauss = By_gauss_N * 0.5 * ( 1. + COS( 3.141592654 * (x-x_0) / (x_max - x_0) ) )
   else
      By_gauss = By_gauss_N * (x/x_0)**3
   end if

END FUNCTION By_gauss
