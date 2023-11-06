!==================================================================================
! This subroutine
!
SUBROUTINE REARRANGE_ACTIVE_BLOCKS

   USE ParallelOperationValues
   USE CurrentProblemValues
   IMPLICIT NONE

   INTEGER s            ! species type

   INTEGER max_buf_len  ! length of bufer, necessary to store the values during the extension of arrays
   INTEGER N_max        ! Maximal number of particles in block after completion of exchange
   INTEGER ALLOC_ERR, DEALLOC_ERR

   REAL(8), ALLOCATABLE :: dbufer(:)  ! bufer for real(8) values
   INTEGER, ALLOCATABLE :: ibufer(:)  ! bufer for integer values

!  INTEGER count        ! counter (cycle) of injected particles

   INTEGER N_start, k

   TYPE(injected_particle), POINTER :: current

   max_buf_len = 0
! find the buffer length if it is necessary to increase the size of some blocks
   DO s = 1, N_spec
      N_max = N_part(s) + N_inject(s)     ! to be sure that we have space for injected particles
      IF (Length(s).GE.N_max) CYCLE
      IF (N_part(s).GT.max_buf_len) max_buf_len = N_part(s)
   END DO
! if necessary, allocate double precision bufer
   ALLOCATE(dbufer(1:max_buf_len), STAT=ALLOC_ERR)
   IF(ALLOC_ERR.NE.0)THEN
      PRINT '(/2x,"Process ",i3," : REARRANGE_ACTIVE_BLOCKS : Error in ALLOCATE dbufer !!!")', Rank_of_process
      PRINT '(2x,"Increasing block size...")'
      PRINT '(2x,"Program will be terminated now :(")'
      STOP
   END IF
! if necessary, allocate integer bufer
   ALLOCATE(ibufer(1:max_buf_len), STAT=ALLOC_ERR)
   IF(ALLOC_ERR.NE.0)THEN
      PRINT '(/2x,"Process ",i3," : REARRANGE_ACTIVE_BLOCKS : Error in ALLOCATE ibufer !!!")', Rank_of_process
      PRINT '(2x,"Increasing block size...")'
      PRINT '(2x,"Program will be terminated now :(")'
      STOP
   END IF
! Increasing the array lengthes (if necessary)
   DO s = 1, N_spec
      N_max = N_part(s) + N_inject(s)

      IF (Length(s).GE.N_max) CYCLE                                    ! if current length is sufficient, go to the other block
      Length(s) = 1.1 * REAL(N_max) + 1                                ! calculate the new length
      IF (N_part(s).EQ.0) THEN                                         ! if block has no particles, just extend the empty block
! X-coordinate
         DEALLOCATE(X_of_spec(s)%part, STAT=DEALLOC_ERR)
         IF(DEALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : Error in DEALLOCATE X_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
            PRINT '(2x,"Increasing block size... No particles...")'
            PRINT '(2x,"Program will be terminated now :(")'
            STOP
         END IF
         ALLOCATE(X_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
         IF(ALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : Error in ALLOCATE X_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
            PRINT '(2x,"Increasing block size... No particles...")'
            PRINT '(2x,"Program will be terminated now :(")'
            STOP
         END IF
         X_of_spec(s)%part(1:Length(s)) = 0.0_8
! X-velocity
         DEALLOCATE(VX_of_spec(s)%part, STAT=DEALLOC_ERR)
         IF(DEALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : Error in DEALLOCATE VX_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
            PRINT '(2x,"Increasing block size... No particles...")'
            PRINT '(2x,"Program will be terminated now :(")'
            STOP
         END IF
         ALLOCATE(VX_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
         IF(ALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : Error in ALLOCATE VX_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
            PRINT '(2x,"Increasing block size... No particles...")'
            PRINT '(2x,"Program will be terminated now :(")'
            STOP
         END IF
         VX_of_spec(s)%part(1:Length(s)) = 0.0_8
! Y-velocity
         DEALLOCATE(VY_of_spec(s)%part, STAT=DEALLOC_ERR)
         IF(DEALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : Error in DEALLOCATE VY_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
            PRINT '(2x,"Increasing block size... No particles...")'
            PRINT '(2x,"Program will be terminated now :(")'
            STOP
         END IF
         ALLOCATE(VY_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
         IF(ALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : Error in ALLOCATE VY_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
            PRINT '(2x,"Increasing block size... No particles...")'
            PRINT '(2x,"Program will be terminated now :(")'
            STOP
         END IF
         VY_of_spec(s)%part(1:Length(s)) = 0.0_8
! Z-velocity
         DEALLOCATE(VZ_of_spec(s)%part, STAT=DEALLOC_ERR)
         IF(DEALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : Error in DEALLOCATE VZ_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
            PRINT '(2x,"Increasing block size... No particles...")'
            PRINT '(2x,"Program will be terminated now :(")'
            STOP
         END IF
         ALLOCATE(VZ_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
         IF(ALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : Error in ALLOCATE VZ_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
            PRINT '(2x,"Increasing block size... No particles...")'
            PRINT '(2x,"Program will be terminated now :(")'
            STOP
         END IF
         VZ_of_spec(s)%part(1:Length(s)) = 0.0_8
! AX-acceleration
         DEALLOCATE(AX_of_spec(s)%part, STAT=DEALLOC_ERR)
         IF(DEALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : Error in DEALLOCATE AX_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
            PRINT '(2x,"Increasing block size... No particles...")'
            PRINT '(2x,"Program will be terminated now :(")'
            STOP
         END IF
         ALLOCATE(AX_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
         IF(ALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : Error in ALLOCATE AX_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
            PRINT '(2x,"Increasing block size... No particles...")'
            PRINT '(2x,"Program will be terminated now :(")'
            STOP
         END IF
         AX_of_spec(s)%part(1:Length(s)) = 0.0_8
! tag
         DEALLOCATE(Tag_of_spec(s)%part, STAT=DEALLOC_ERR)
         IF(DEALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : Error in DEALLOCATE Tag_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
            PRINT '(2x,"Increasing block size... No particles...")'
            PRINT '(2x,"Program will be terminated now :(")'
            STOP
         END IF
         ALLOCATE(Tag_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
         IF(ALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : Error in ALLOCATE Tag_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
            PRINT '(2x,"Increasing block size... No particles...")'
            PRINT '(2x,"Program will be terminated now :(")'
            STOP
         END IF
         Tag_of_spec(s)%part(1:Length(s)) = 0
! prev-X-velocity
         DEALLOCATE(prev_VX_of_spec(s)%part, STAT=DEALLOC_ERR)
         IF(DEALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : Error in DEALLOCATE prev_VX_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
            PRINT '(2x,"Increasing block size... No particles...")'
            PRINT '(2x,"Program will be terminated now :(")'
            STOP
         END IF
         ALLOCATE(prev_VX_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
         IF(ALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : Error in ALLOCATE prev_VX_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
            PRINT '(2x,"Increasing block size... No particles...")'
            PRINT '(2x,"Program will be terminated now :(")'
            STOP
         END IF
         prev_VX_of_spec(s)%part(1:Length(s)) = 0.0_8

         CYCLE        ! take the next species

      END IF

! if we are here, then the block has particles and it must be extended, therefore we have to store the data at
! the bufer at first, then extend the array (deallocate-allocate), and then copy the data back from the bufer

! X-coordinate, with bufer
      dbufer(1:N_part(s)) = X_of_spec(s)%part(1:N_part(s))
      DEALLOCATE(X_of_spec(s)%part, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : Error in DEALLOCATE X_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
         PRINT '(2x,"Increasing block size...")'
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF
      ALLOCATE(X_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : Error in ALLOCATE X_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
         PRINT '(2x,"Increasing block size...")'
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF
      X_of_spec(s)%part(1:N_part(s)) = dbufer(1:N_part(s))
      X_of_spec(s)%part((N_part(s)+1):Length(s)) = 0.0_8
! X-velocity, with bufer
      dbufer(1:N_part(s)) = VX_of_spec(s)%part(1:N_part(s))
      DEALLOCATE(VX_of_spec(s)%part, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : Error in DEALLOCATE VX_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
         PRINT '(2x,"Increasing block size...")'
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF
      ALLOCATE(VX_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : Error in ALLOCATE VX_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
         PRINT '(2x,"Increasing block size...")'
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF
      VX_of_spec(s)%part(1:N_part(s)) = dbufer(1:N_part(s))
      VX_of_spec(s)%part((N_part(s)+1):Length(s)) = 0.0_8
! Y-velocity, with bufer
      dbufer(1:N_part(s)) = VY_of_spec(s)%part(1:N_part(s))
      DEALLOCATE(VY_of_spec(s)%part, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : Error in DEALLOCATE VY_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
         PRINT '(2x,"Increasing block size...")'
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF
      ALLOCATE(VY_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : Error in ALLOCATE VY_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
         PRINT '(2x,"Increasing block size...")'
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF
      VY_of_spec(s)%part(1:N_part(s)) = dbufer(1:N_part(s))
      VY_of_spec(s)%part((N_part(s)+1):Length(s)) = 0.0_8
! Z-velocity, with bufer
      dbufer(1:N_part(s)) = VZ_of_spec(s)%part(1:N_part(s))
      DEALLOCATE(VZ_of_spec(s)%part, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : Error in DEALLOCATE VZ_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
         PRINT '(2x,"Increasing block size...")'
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF
      ALLOCATE(VZ_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : Error in ALLOCATE VZ_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
         PRINT '(2x,"Increasing block size...")'
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF
      VZ_of_spec(s)%part(1:N_part(s)) = dbufer(1:N_part(s))
      VZ_of_spec(s)%part((N_part(s)+1):Length(s)) = 0.0_8
! AX-acceleration, with bufer
      dbufer(1:N_part(s)) = AX_of_spec(s)%part(1:N_part(s))
      DEALLOCATE(AX_of_spec(s)%part, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : Error in DEALLOCATE AX_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
         PRINT '(2x,"Increasing block size...")'
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF
      ALLOCATE(AX_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : Error in ALLOCATE AX_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
         PRINT '(2x,"Increasing block size...")'
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF
      AX_of_spec(s)%part(1:N_part(s)) = dbufer(1:N_part(s))
      AX_of_spec(s)%part((N_part(s)+1):Length(s)) = 0.0_8
! tag, with bufer
      ibufer(1:N_part(s)) = Tag_of_spec(s)%part(1:N_part(s))
      DEALLOCATE(Tag_of_spec(s)%part, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : Error in DEALLOCATE Tag_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
         PRINT '(2x,"Increasing block size...")'
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF
      ALLOCATE(Tag_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : Error in ALLOCATE Tag_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
         PRINT '(2x,"Increasing block size...")'
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF
      Tag_of_spec(s)%part(1:N_part(s)) = ibufer(1:N_part(s))
      Tag_of_spec(s)%part((N_part(s)+1):Length(s)) = 0
! prev-X-velocity, with bufer
      dbufer(1:N_part(s)) = prev_VX_of_spec(s)%part(1:N_part(s))
      DEALLOCATE(prev_VX_of_spec(s)%part, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : Error in DEALLOCATE prev_VX_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
         PRINT '(2x,"Increasing block size...")'
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF
      ALLOCATE(prev_VX_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : Error in ALLOCATE prev_VX_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
         PRINT '(2x,"Increasing block size...")'
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF
      prev_VX_of_spec(s)%part(1:N_part(s)) = dbufer(1:N_part(s))
      prev_VX_of_spec(s)%part((N_part(s)+1):Length(s)) = 0.0_8
   END DO

! deallocate double precision bufer
   IF (ALLOCATED(dbufer)) THEN
      DEALLOCATE(dbufer, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(/2x,"Process ",i3," : REARRANGE_ACTIVE_BLOCKS : Error in DEALLOCATE dbufer !!!")', Rank_of_process
         PRINT '(2x,"Increasing block size...")'
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF
   END IF
! deallocate integer bufer
   IF(ALLOCATED(ibufer)) THEN
      DEALLOCATE(ibufer, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(/2x,"Process ",i3," : REARRANGE_ACTIVE_BLOCKS : Error in DEALLOCATE ibufer !!!")', Rank_of_process
         PRINT '(2x,"Increasing block size...")'
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF
   END IF
! end of arrays extension
!===========================
!print *, 'started recursive exchange', T_cntr
!print *, N_accept
!print *, Counter_add

   DO s = 1, N_spec

      NULLIFY(current)

      SELECT CASE(s)
       CASE(1)                                 ! inject electrons into the first group
         current => Inj_electron
       CASE(2)                                 ! inject ions to the first group
         current => Inj_ion
      END SELECT

      N_start = N_part(s) + 1
      N_part(s) = N_part(s) + N_inject(s)

      DO k = N_start, N_part(s)

         IF(.NOT.ASSOCIATED(current)) THEN
            PRINT '(2x,"Process ",i3," : REARRANGE_ACTIVE_BLOCKS : While adding injected species ",i2," ...")', Rank_of_process, s
            PRINT '(2x,"Error! current is NOT associated!")'
            PRINT '(2x,"Program will be terminated now, sorry! :(")'
            STOP
         END IF
         X_of_spec(s)%part(k) = current%X
         VX_of_spec(s)%part(k) = current%VX
         VY_of_spec(s)%part(k) = current%VY
         VZ_of_spec(s)%part(k) = current%VZ
         AX_of_spec(s)%part(k) = current%AX
         Tag_of_spec(s)%part(k) = current%Tag
         IF (k.LT.N_part(s)) current => current%next
      END DO

      IF (N_part(s).GT.Length(s)) THEN
         PRINT '(2x,"Process ",i3," : REARRANGE_ACTIVE_BLOCKS : While rearranging species ",i2," ...")', Rank_of_process, s
         PRINT '(2x,"Error! Number of particles is greater than the length!")'
         PRINT '(2x," N_part = ",i7," Length = ",i7)', N_part(s), Length(s)
         PRINT '(2x,"Program will be terminated now, sorry! :(")'
         STOP
      END IF

   END DO

!print *, 'finished recursive exchange'
!===========================

   max_buf_len = 0
! find the buffer length if it is necessary to decrease the size of some blocks
   DO s = 1, N_spec
      N_max = 0.8*REAL(Length(s))
      IF (N_part(s).GE.N_max) CYCLE
      IF (N_part(s).GT.max_buf_len) max_buf_len = N_part(s)   ! if for some species the particles occupy less than 80% of the assigned array
   END DO

! if necessary, allocate double precision bufer
   ALLOCATE(dbufer(1:max_buf_len), STAT=ALLOC_ERR)
   IF(ALLOC_ERR.NE.0)THEN
      PRINT '(/2x,"Process ",i3," : REARRANGE_ACTIVE_BLOCKS : Error in ALLOCATE dbufer !!!")', Rank_of_process
      PRINT '(2x,"Decreasing block size...")'
      PRINT '(2x,"Program will be terminated now :(")'
      STOP
   END IF
! if necessary, allocate integer bufer
   ALLOCATE(ibufer(1:max_buf_len), STAT=ALLOC_ERR)
   IF(ALLOC_ERR.NE.0)THEN
      PRINT '(/2x,"Process ",i3," : REARRANGE_ACTIVE_BLOCKS : Error in ALLOCATE ibufer !!!")', Rank_of_process
      PRINT '(2x,"Decreasing block size...")'
      PRINT '(2x,"Program will be terminated now :(")'
      STOP
   END IF

! Decreasing the array lengthes, if necessary
   DO s = 1, N_spec
      N_max = 0.8*REAL(Length(s))
      IF (N_part(s).GE.N_max) CYCLE     ! if block is not too lengthy, take the next block
      IF (Length(s).LT.10) CYCLE     ! to avoid cycling 10-9-10 or 9-8-9, etc.

      Length(s) = 0.9*REAL(Length(s)) + 1   ! this is the new length

      IF (N_part(s).GT.Length(s)) THEN
         PRINT '(2x,"Process ",i3," : While decreasing arrays...")', Rank_of_process
         PRINT '(2x,"Error! Number of particles is greater than the length!")'
         PRINT '(2x,"Species #",i5," N_part = ",i7," Length = ",i7)', s, N_part(s), Length(s)
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF

      IF (N_part(s).GT.0) THEN          ! if the block has particles
! X-coordinate
         dbufer(1:N_part(s)) = X_of_spec(s)%part(1:N_part(s))
         DEALLOCATE(X_of_spec(s)%part, STAT=DEALLOC_ERR)
         IF(DEALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : Error in DEALLOCATE X_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
            PRINT '(2x,"Decreasing block size...")'
            PRINT '(2x,"Program will be terminated now :(")'
            STOP
         END IF
         ALLOCATE(X_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
         IF(ALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : Error in ALLOCATE X_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
            PRINT '(2x,"Decreasing block size...")'
            PRINT '(2x,"Program will be terminated now :(")'
            STOP
         END IF
         X_of_spec(s)%part(1:N_part(s)) = dbufer(1:N_part(s))
         X_of_spec(s)%part((N_part(s)+1):Length(s)) = 0.0_8
! X-velocity
         dbufer(1:N_part(s)) = VX_of_spec(s)%part(1:N_part(s))
         DEALLOCATE(VX_of_spec(s)%part, STAT=DEALLOC_ERR)
         IF(DEALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : Error in DEALLOCATE VX_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
            PRINT '(2x,"Decreasing block size...")'
            PRINT '(2x,"Program will be terminated now :(")'
            STOP
         END IF
         ALLOCATE(VX_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
         IF(ALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : Error in ALLOCATE VX_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
            PRINT '(2x,"Increasing block size...")'
            PRINT '(2x,"Program will be terminated now :(")'
            STOP
         END IF
         VX_of_spec(s)%part(1:N_part(s)) = dbufer(1:N_part(s))
         VX_of_spec(s)%part((N_part(s)+1):Length(s)) = 0.0_8
! Y-velocity
         dbufer(1:N_part(s)) = VY_of_spec(s)%part(1:N_part(s))
         DEALLOCATE(VY_of_spec(s)%part, STAT=DEALLOC_ERR)
         IF(DEALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : Error in DEALLOCATE VY_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
            PRINT '(2x,"Decreasing block size...")'
            PRINT '(2x,"Program will be terminated now :(")'
            STOP
         END IF
         ALLOCATE(VY_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
         IF(ALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : Error in ALLOCATE VY_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
            PRINT '(2x,"Decreasing block size...")'
            PRINT '(2x,"Program will be terminated now :(")'
            STOP
         END IF
         VY_of_spec(s)%part(1:N_part(s)) = dbufer(1:N_part(s))
         VY_of_spec(s)%part((N_part(s)+1):Length(s)) = 0.0_8
! Z-velocity
         dbufer(1:N_part(s)) = VZ_of_spec(s)%part(1:N_part(s))
         DEALLOCATE(VZ_of_spec(s)%part, STAT=DEALLOC_ERR)
         IF(DEALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : Error in DEALLOCATE VZ_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
            PRINT '(2x,"Decreasing block size...")'
            PRINT '(2x,"Program will be terminated now :(")'
            STOP
         END IF
         ALLOCATE(VZ_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
         IF(ALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : Error in ALLOCATE VZ_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
            PRINT '(2x,"Decreasing block size...")'
            PRINT '(2x,"Program will be terminated now :(")'
            STOP
         END IF
         VZ_of_spec(s)%part(1:N_part(s)) = dbufer(1:N_part(s))
         VZ_of_spec(s)%part((N_part(s)+1):Length(s)) = 0.0_8
! AX-acceleration
         dbufer(1:N_part(s)) = AX_of_spec(s)%part(1:N_part(s))
         DEALLOCATE(AX_of_spec(s)%part, STAT=DEALLOC_ERR)
         IF(DEALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : Error in DEALLOCATE AX_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
            PRINT '(2x,"Decreasing block size...")'
            PRINT '(2x,"Program will be terminated now :(")'
            STOP
         END IF
         ALLOCATE(AX_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
         IF(ALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : Error in ALLOCATE AX_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
            PRINT '(2x,"Decreasing block size...")'
            PRINT '(2x,"Program will be terminated now :(")'
            STOP
         END IF
         AX_of_spec(s)%part(1:N_part(s)) = dbufer(1:N_part(s))
         AX_of_spec(s)%part((N_part(s)+1):Length(s)) = 0.0_8
! tag
         ibufer(1:N_part(s)) = Tag_of_spec(s)%part(1:N_part(s))
         DEALLOCATE(Tag_of_spec(s)%part, STAT=DEALLOC_ERR)
         IF(DEALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : Error in DEALLOCATE Tag_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
            PRINT '(2x,"Decreasing block size...")'
            PRINT '(2x,"Program will be terminated now :(")'
            STOP
         END IF
         ALLOCATE(Tag_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
         IF(ALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : Error in ALLOCATE Tag_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
            PRINT '(2x,"Decreasing block size...")'
            PRINT '(2x,"Program will be terminated now :(")'
            STOP
         END IF
         Tag_of_spec(s)%part(1:N_part(s)) = ibufer(1:N_part(s))
         Tag_of_spec(s)%part((N_part(s)+1):Length(s)) = 0
! prev-X-velocity
         dbufer(1:N_part(s)) = prev_VX_of_spec(s)%part(1:N_part(s))
         DEALLOCATE(prev_VX_of_spec(s)%part, STAT=DEALLOC_ERR)
         IF(DEALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : Error in DEALLOCATE prev_VX_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
            PRINT '(2x,"Decreasing block size...")'
            PRINT '(2x,"Program will be terminated now :(")'
            STOP
         END IF
         ALLOCATE(prev_VX_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
         IF(ALLOC_ERR.NE.0)THEN
            PRINT '(2x,"Process ",i3," : Error in ALLOCATE prev_VX_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
            PRINT '(2x,"Increasing block size...")'
            PRINT '(2x,"Program will be terminated now :(")'
            STOP
         END IF
         prev_VX_of_spec(s)%part(1:N_part(s)) = dbufer(1:N_part(s))
         prev_VX_of_spec(s)%part((N_part(s)+1):Length(s)) = 0.0_8

         CYCLE                                                      ! take the next species

      END IF

! if we are here, then the number of particles in block is equal to zero, but the length is still large
! therefore we shall just deallocate/allocate the clean arrays

! X-coordinate
      DEALLOCATE(X_of_spec(s)%part, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : Error in DEALLOCATE X_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
         PRINT '(2x,"Decreasing block size... No particles...")'
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF
      ALLOCATE(X_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : Error in ALLOCATE X_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
         PRINT '(2x,"Decreasing block size... No particles...")'
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF
      X_of_spec(s)%part(1:Length(s)) = 0.0_8
! X-velocity
      DEALLOCATE(VX_of_spec(s)%part, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : Error in DEALLOCATE VX_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
         PRINT '(2x,"Decreasing block size... No particles...")'
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF
      ALLOCATE(VX_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : Error in ALLOCATE VX_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
         PRINT '(2x,"Decreasing block size... No particles...")'
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF
      VX_of_spec(s)%part(1:Length(s)) = 0.0_8
! Y-velocity
      DEALLOCATE(VY_of_spec(s)%part, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : Error in DEALLOCATE VY_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
         PRINT '(2x,"Decreasing block size... No particles...")'
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF
      ALLOCATE(VY_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : Error in ALLOCATE VY_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
         PRINT '(2x,"Decreasing block size... No particles...")'
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF
      VY_of_spec(s)%part(1:Length(s)) = 0.0_8
! Z-velocity
      DEALLOCATE(VZ_of_spec(s)%part, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : Error in DEALLOCATE VZ_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
         PRINT '(2x,"Decreasing block size... No particles...")'
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF
      ALLOCATE(VZ_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : Error in ALLOCATE VZ_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
         PRINT '(2x,"Decreasing block size... No particles...")'
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF
      VZ_of_spec(s)%part(1:Length(s)) = 0.0_8
! AX-acceleration
      DEALLOCATE(AX_of_spec(s)%part, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : Error in DEALLOCATE AX_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
         PRINT '(2x,"Decreasing block size... No particles...")'
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF
      ALLOCATE(AX_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : Error in ALLOCATE AX_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
         PRINT '(2x,"Decreasing block size... No particles...")'
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF
      AX_of_spec(s)%part(1:Length(s)) = 0.0_8
! tag
      DEALLOCATE(Tag_of_spec(s)%part, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : Error in DEALLOCATE Tag_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
         PRINT '(2x,"Decreasing block size... No particles...")'
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF
      ALLOCATE(Tag_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : Error in ALLOCATE Tag_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
         PRINT '(2x,"Decreasing block size... No particles...")'
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF
      Tag_of_spec(s)%part(1:Length(s)) = 0
! prev-X-velocity
      DEALLOCATE(prev_VX_of_spec(s)%part, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : Error in DEALLOCATE prev_VX_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
         PRINT '(2x,"Decreasing block size... No particles...")'
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF
      ALLOCATE(prev_VX_of_spec(s)%part(1:Length(s)), STAT=ALLOC_ERR)
      IF(ALLOC_ERR.NE.0)THEN
         PRINT '(2x,"Process ",i3," : Error in ALLOCATE prev_VX_of_spec(s)%part !!! s = ",i2)', Rank_of_process, s
         PRINT '(2x,"Decreasing block size... No particles...")'
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF
      prev_VX_of_spec(s)%part(1:Length(s)) = 0.0_8
   END DO

! deallocate double precision bufer
   IF (ALLOCATED(dbufer)) THEN
      DEALLOCATE(dbufer, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(/2x,"Process ",i3," : REARRANGE_ACTIVE_BLOCKS : Error in DEALLOCATE dbufer !!!")', Rank_of_process
         PRINT '(2x,"Decreasing block size...")'
         PRINT  '(2x,"Program will be terminated now :(")'
         STOP
      END IF
   END IF
! deallocate integer bufer
   IF (ALLOCATED(ibufer)) THEN
      DEALLOCATE(ibufer, STAT=DEALLOC_ERR)
      IF(DEALLOC_ERR.NE.0)THEN
         PRINT '(/2x,"Process ",i3," : REARRANGE_ACTIVE_BLOCKS : Error in DEALLOCATE ibufer !!!")', Rank_of_process
         PRINT '(2x,"Decreasing block size...")'
         PRINT '(2x,"Program will be terminated now :(")'
         STOP
      END IF
   END IF
! end of arrays extension

!print *, 'finished REARRANGE_ACTIVE_BLOCK'

END SUBROUTINE  REARRANGE_ACTIVE_BLOCKS
