!===================================================================================================
SUBROUTINE INITIATE_MC_COLLISIONS

  use mpi
  USE ParallelOperationValues
  USE MCCollisions
  USE CurrentProblemValues, ONLY : N_max_vel, T_e_eV, T_i_eV, N_spec, Ms, N_box_vel, r_max_vel

  IMPLICIT NONE

!  INCLUDE 'mpif.h'

  LOGICAL exists
  INTEGER s, k

  CHARACTER (67) buf
  INTEGER ierr

  maxNcolkind_spec(1) = 5            ! Currently we have 5 kinds of e-n collisions         !@#$
  maxNcolkind_spec(2) = 3            ! Currently we have 3 kinds of i-n collisions         !@#$
! maxNcolkind_spec(3) = 0
  
  Colflag_kind_spec = 0
 
  INQUIRE (FILE = 'ssc_partcolls.dat', EXIST = exists)
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  OPEN (9, FILE = 'ssc_partcolls.dat')

  IF(exists) THEN

     PRINT '(2x,"Process ",i3," : ssc_partcolls.dat file is found. Reading the data file...")', Rank_of_process

     READ (9, '(A67)') buf !================= NEUTRAL COMPONENT PARAMETERS ====================")')
     READ (9, '(A67)') buf !-------d----- Neutral gas species ( 1 = Argon, 0 = Helium ) -------")')
     READ (9, '(7x,i1)') Neutral_flag
     READ (9, '(A67)') buf !--dddddd.ddd- Mass (a.m.u.) ---------------------------------------")')
     READ (9, '(2x,f10.3)') M_neutral_amu 
     READ (9, '(A67)') buf !--#d.dddE#dd- Density (m^-3) --------------------------------------")')
     READ (9, '(2x,e10.3)') N_neutral_m3 
     READ (9, '(A67)') buf !--dddddd.ddd- Temperature (eV) ------------------------------------")')
     READ (9, '(2x,f10.3)') T_neutral_eV     
     READ (9, '(A67)') buf !============ ELECTRON - NEUTRAL COLLISIONS, ACTIVATION ============")')
     READ (9, '(A67)') buf !-------d----- Elastic-1 (1 = yes, 0 = no) -------------------------")')
     READ (9, '(7x,i1)') Colflag_kind_spec(1,1)
     READ (9, '(A67)') buf !-------d----- Excitation-1 (1 = yes, 0 = no) ----------------------")')
     READ (9, '(7x,i1)') Colflag_kind_spec(2,1)
     READ (9, '(A67)') buf !-------d----- Excitation-2 (1 = yes, 0 = no) ----------------------")')
     READ (9, '(7x,i1)') Colflag_kind_spec(5,1)
     READ (9, '(A67)') buf !-------d----- Ionization-1 (1 = yes, 0 = no) ----------------------")')
     READ (9, '(7x,i1)') Colflag_kind_spec(3,1)
     READ (9, '(A67)') buf !============== ELECTRON - NEUTRAL COLLISIONS, CONTROL =============")')
     READ (9, '(A67)') buf !--dddddd.ddd- Maximal electron energy (eV), default if negative ---")')
     READ (9, '(2x,f10.3)') maxEnergy_eV_spec(1) 
     READ (9, '(A67)') buf !--dddddd----- Number of energy values (>0) ------------------------")')
     READ (9, '(2x,i6)') Nenergyval_spec(1)
     READ (9, '(A67)') buf !========== ELECTRON - TURBULENCE, ACTIVATION and CONTROL ==========")')         !@#$
     READ (9, '(A67)') buf !-------d----- Turbulence-1 (1 = yes, 0 = no) ----------------------")')         !@#$
     READ (9, '(7x,i1)') Colflag_kind_spec(4,1)                                                            !@#$
     READ (9, '(A67)') buf !--#d.dddE#dd- Frequency, model-1 (s^-1) ---------------------------")')         !@#$      
     READ (9, '(2x,e10.3)') Freq_turb_e_1_s1                                                                  !@#$
     READ (9, '(A67)') buf !============== ION - NEUTRAL COLLISIONS, ACTIVATION ===============")')
     READ (9, '(A67)') buf !-------d----- Elastic-1 (1 = yes, 0 = no) -------------------------")')
     READ (9, '(7x,i1)') Colflag_kind_spec(1,2)
     READ (9, '(A67)') buf !-------d----- Charge exchange-1 (1 = yes, 0 = no) -----------------")')
     READ (9, '(7x,i1)') Colflag_kind_spec(2,2)
     READ (9, '(A67)') buf !================ ION - NEUTRAL COLLISIONS, CONTROL ================")')
     READ (9, '(A67)') buf !--dddddd.ddd- Maximal ion energy (eV), default if negative --------")')
     READ (9, '(2x,f10.3)') maxEnergy_eV_spec(2)
     READ (9, '(A67)') buf !--dddddd----- Number of energy values (>0) ------------------------")')
     READ (9, '(2x,i6)') Nenergyval_spec(2)
     READ (9, '(A67)') buf !============ ION - TURBULENCE, ACTIVATION and CONTROL =============")')         !@#$
     READ (9, '(A67)') buf !-------d----- Turbulence-1 (1 = yes, 0 = no) ----------------------")')         !@#$
     READ (9, '(7x,i1)') Colflag_kind_spec(3,2)                                                            !@#$
     READ (9, '(A67)') buf !--#d.dddE#dd- Frequency, model-1 (s^-1) ---------------------------")')         !@#$      
     READ (9, '(2x,e10.3)') Freq_turb_i_1_s1                                                                  !@#$
!     Colflag_kind_spec(4,2) = 1
     
  ELSE

     Neutral_flag           = 1           ! 1 = Argon, 0 = Helium
     M_neutral_amu          = 131.0_8
     N_neutral_m3           = 2.0e18
     T_neutral_eV           = 0.1_8 
     Colflag_kind_spec(1,1) = 1           ! e-n, elastic-1,       yes
     Colflag_kind_spec(2,1) = 1           ! e-n, excitation-1,    yes
     Colflag_kind_spec(5,1) = 0           ! e-n, excitation-2,    yes
     Colflag_kind_spec(3,1) = 1           ! e-n, ionization-1,    yes
     maxEnergy_eV_spec(1)   = 1000.0_8    ! 
     Nenergyval_spec(1)     = 20001
     Colflag_kind_spec(4,1) = 1           ! e-turbulence-1,       no                                !@#$
     Freq_turb_e_1_s1       = -80.0_8                                                                 !@#$
     Colflag_kind_spec(1,2) = 0           ! i-n, elastic-1,         no
     Colflag_kind_spec(2,2) = 0           ! i-n, charge exchange-1, no
     maxEnergy_eV_spec(2)   = -1.0_8      ! default value will be calculated
     Nenergyval_spec(2)     = 101
     Colflag_kind_spec(3,2) = 0           ! i-turbulence-1,       no                                !@#$
     Freq_turb_i_1_s1       = 0.0_8                                                                 !@#$
     Colflag_kind_spec(4,2) = 0           ! empty                                                   !@#$

     PRINT '(2x,"Process ",i3," : File with the name ssc_partcolls.dat not found. Use the default settings ...")', &
                                                                                                    & Rank_of_process

     IF (Rank_of_process.EQ.0) THEN

        PRINT '(/2x,"Process ",i3," : Create ssc_partcolls.dat file . . .")', Rank_of_process

        WRITE (9, '("================= NEUTRAL COMPONENT PARAMETERS ====================")')
        WRITE (9, '("-------d----- Neutral gas species ( 1 = Argon, 0 = Helium ) -------")')
        WRITE (9, '(7x,i1)') Neutral_flag
        WRITE (9, '("--dddddd.ddd- Mass (a.m.u.) ---------------------------------------")')
        WRITE (9, '(2x,f10.3)') M_neutral_amu 
        WRITE (9, '("--#d.dddE#dd- Density (m^-3) --------------------------------------")')
        WRITE (9, '(2x,e10.3)') N_neutral_m3 
        WRITE (9, '("--dddddd.ddd- Temperature (eV) ------------------------------------")')
        WRITE (9, '(2x,f10.3)') T_neutral_eV     
        WRITE (9, '("============ ELECTRON - NEUTRAL COLLISIONS, ACTIVATION ============")')
        WRITE (9, '("-------d----- Elastic-1 (1 = yes, 0 = no) -------------------------")')
        WRITE (9, '(7x,i1)') Colflag_kind_spec(1,1)
        WRITE (9, '("-------d----- Excitation-1 (1 = yes, 0 = no) ----------------------")')
        WRITE (9, '(7x,i1)') Colflag_kind_spec(2,1)
        WRITE (9, '("-------d----- Excitation-1 (1 = yes, 0 = no) ----------------------")')
        WRITE (9, '(7x,i1)') Colflag_kind_spec(5,1)
        WRITE (9, '("-------d----- Ionization-1 (1 = yes, 0 = no) ----------------------")')
        WRITE (9, '(7x,i1)') Colflag_kind_spec(3,1)
        WRITE (9, '("============== ELECTRON - NEUTRAL COLLISIONS, CONTROL =============")')
        WRITE (9, '("--dddddd.ddd- Maximal electron energy (eV), default if negative ---")')
        WRITE (9, '(2x,f10.3)') maxEnergy_eV_spec(1) 
        WRITE (9, '("--dddddd----- Number of energy values (>0) ------------------------")')
        WRITE (9, '(2x,i6)') Nenergyval_spec(1)
        WRITE (9, '("========== ELECTRON - TURBULENCE, ACTIVATION and CONTROL ==========")')         !@#$
        WRITE (9, '("-------d----- Turbulence-1 (1 = yes, 0 = no) ----------------------")')         !@#$
        WRITE (9, '(7x,i1)') Colflag_kind_spec(4,1)                                                            !@#$
        WRITE (9, '("--#d.dddE#dd- Frequency, model-1 (s^-1) ---------------------------")')         !@#$      
        WRITE (9, '(2x,e10.3)') Freq_turb_e_1_s1                                                                  !@#$
        WRITE (9, '("============== ION - NEUTRAL COLLISIONS, ACTIVATION ===============")')
        WRITE (9, '("-------d----- Elastic-1 (1 = yes, 0 = no) -------------------------")')
        WRITE (9, '(7x,i1)') Colflag_kind_spec(1,2)
        WRITE (9, '("-------d----- Charge exchange-1 (1 = yes, 0 = no) -----------------")')
        WRITE (9, '(7x,i1)') Colflag_kind_spec(2,2)
        WRITE (9, '("================ ION - NEUTRAL COLLISIONS, CONTROL ================")')
        WRITE (9, '("--dddddd.ddd- Maximal ion energy (eV), default if negative --------")')
        WRITE (9, '(2x,f10.3)') maxEnergy_eV_spec(2)
        WRITE (9, '("--dddddd----- Number of energy values (>0) ------------------------")')
        WRITE (9, '(2x,i6)') Nenergyval_spec(2)
        WRITE (9, '("============ ION - TURBULENCE, ACTIVATION and CONTROL =============")')         !@#$
        WRITE (9, '("-------d----- Turbulence-1 (1 = yes, 0 = no) ----------------------")')         !@#$
        WRITE (9, '(7x,i1)') Colflag_kind_spec(3,2)                                                            !@#$
        WRITE (9, '("--#d.dddE#dd- Frequency, model-1 (s^-1) ---------------------------")')         !@#$      
        WRITE (9, '(2x,e10.3)') Freq_turb_i_1_s1                                                                  !@#$
     END IF

  END IF

  CLOSE (9, STATUS = 'KEEP')

  IF (maxEnergy_eV_spec(1).LT.0.0_8) THEN                                 ! 
     maxEnergy_eV_spec(1) = 2.0_8 * r_max_vel * r_max_vel * T_e_eV        ! default value of the maximal energy for electrons
  END IF                                                                  !

  IF (maxEnergy_eV_spec(2).LT.0.0_8) THEN                                 ! 
     maxEnergy_eV_spec(2) = 2.0_8 * r_max_vel * r_max_vel * T_i_eV        ! default value of the maximal energy for ions
  END IF                                                                  !
  
!  IF (Freq_turb_e_1_s1.LT.0.0_8) THEN                                     ! frequency of electron turbulence collisions is defined
!     Freq_turb_e_1_s1 = ABS(W_cycl_x_s1 / Freq_turb_e_1_s1)                 ! by the electron cyclotron frequency
!  END IF                                                                  ! (if the input data is negative)

!print *, '1 <- ', Rank_of_process
  CALL CONFIG_READ_CRSECT_ARRAYS
!print *, '2 <- ', Rank_of_process
  CALL CONFIGURE_COLLSPEC_ARRAYS
  CALL SETVALUES_COLLSPEC_ARRAYS
!print *, '3 <- ', Rank_of_process
  CALL CONFIGURE_COLLISION_ARRAYS
  CALL SETVALUES_COLLISION_ARRAYS
!print *, '4 <- ', Rank_of_process

  CALL PrepareMaxwellDistribIntegral

  IF (Rank_of_process.EQ.0) THEN

     DO s = 1, N_spec
        SELECT CASE (s)
           CASE (1)
              PRINT '(/2x,"Activated electron-neutral/turbulence collisions:")'
              IF (Ncolkind_spec(s).EQ.0) THEN
                 PRINT '(4x,"None :)")'
                 CYCLE
              END IF
              DO k = 1, Ncolkind_spec(s)
                 SELECT CASE (Colkind_of_spec(s)%activated(k))
                    CASE(1)  
                       PRINT '(4x,"Elastic, model 1")'
                    CASE(2)
                       PRINT '(4x,"Excitation, model 1")'
                    CASE(3)
                       PRINT '(4x,"Ionization, model 1")'
                    CASE(4)                                                                                                       !@#$
                       PRINT '(4x,"Turbulence, model 1, frequency is ",e10.3," s^-1")', Freq_turb_e_1_s1                          !@#$
                    CASE(5)                                                                                                       !@#$ 
                       PRINT '(4x,"Excitation, model 2")'                                                                         !@#$
                 END SELECT
              END DO

           CASE (2)
              PRINT '(/2x,"Activated ion-neutral/turbulence collisions:")'
              IF (Ncolkind_spec(s).EQ.0) THEN
                 PRINT '(4x,"None :)")'
                 CYCLE
              END IF
              DO k = 1, Ncolkind_spec(s)
                 SELECT CASE (Colkind_of_spec(s)%activated(k))
                    CASE(1)
                       PRINT '(4x,"Elastic, model 1")'
                    CASE(2)
                       PRINT '(4x,"Charge exchange, model 1")'
                    CASE(3)                                                                                                       !@#$
                       PRINT '(4x,"Turbulence, model 1, frequency is ",e10.3," s^-1")', Freq_turb_i_1_s1                          !@#$
                  END SELECT
               END DO

        END SELECT
        PRINT '(/2x,"Colliding portion (%): ",f8.5)', 100.0 * maxColfrac_spec(s)
     END DO

     SELECT CASE (Neutral_flag)
        CASE (0)
          PRINT '(/2x,"Using Helium collision parameters for CollideElectron_# and CollideIon_# functions.")'
          ! Check if the inputted masses are consistent with the neutral flag
          IF ((M_neutral_amu.LT.4.000_8).OR.(M_neutral_amu.GT.4.020_8)) THEN ! boundaries chosen without a specific tolerance in mind
            PRINT '(/2x,"Error! Mass of the neutral gas is not consistent with the neutral flag (helium).")'
            PRINT '(/2x,"You inputted: ",f10.3," a.m.u.")', M_neutral_amu
            PRINT '(2x,"Program will run, but be warned :(")'
          END IF

        CASE (1)
          PRINT '(/2x,"Using Argon collision parameters for CollideElectron_# and CollideIon_# functions.")'
          ! Check if the inputted masses are consistent with the neutral flag
          IF ((M_neutral_amu.LT.39.930_8).OR.(M_neutral_amu.GT.39.960)) THEN ! boundaries chosen without a specific tolerance in mind
            PRINT '(/2x,"Error! Mass of the neutral gas is not consistent with the neutral flag (argon).")'
            PRINT '(/2x,"You inputted: ",f10.3," a.m.u.")', M_neutral_amu
            PRINT '(2x,"Program will run, but be warned :(")'
          END IF

        CASE DEFAULT 
          PRINT '(/2x,"ERROR: No case met in Neutral_flag!!!")'
          STOP

     END SELECT
  END IF

  e_n_1_count = 0; e_n_2_count = 0; e_n_3_count = 0; i_n_1_count = 0; i_n_2_count = 0
  e_t_4_count = 0; i_t_3_count = 0                                                              !@#$
  e_n_5_count = 0
! **** 0.00054858 = m_e_kg / 1_amu_kg = 9.109534e-31 / 1.660565e-27 
  alpha_Vscl = SQRT(0.00054858_8 * T_neutral_eV / (T_e_eV * M_neutral_amu)) 
  
  DO s = 1, N_spec
     energy_factor_eV_s(s) = T_e_eV * Ms(s) / N_box_vel**2
  END DO
   
!  PRINT '(2x,"Finished :)")'

END SUBROUTINE INITIATE_MC_COLLISIONS 

!=========================================================================================================
! Note, this subroutine must be called before any first call of the FREQUENCY_OF_COLLISION function
!
SUBROUTINE CONFIG_READ_CRSECT_ARRAYS 

  USE ParallelOperationValues
  USE MCCollisions

  IMPLICIT NONE

  INTEGER ALLOC_ERR
  LOGICAL exists
  INTEGER j

! read cross sections of e-n elastic collisions from data file, if this kind of collisions is activated   
  IF (Colflag_kind_spec(1, 1).EQ.1) THEN

     INQUIRE (FILE = 'ssc_crsect_en_elast.dat', EXIST = exists)
     
     IF (exists) THEN
     
        OPEN (9, FILE = 'ssc_crsect_en_elast.dat')

        PRINT '(2x,"Process ",i3," : e-n elastic collisions cross-sections data file is found. Reading the data file...")', &
                                                                                                              &Rank_of_process
        READ (9, '(2x,i4)') N_en_elast

        ALLOCATE(Energy_en_elast_eV(1:N_en_elast), STAT=ALLOC_ERR)
        IF(ALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in ALLOCATE Energy_en_elast_eV !!!")', Rank_of_process
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF

        ALLOCATE(CrSect_en_elast_m2(1:N_en_elast), STAT=ALLOC_ERR)
        IF(ALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in ALLOCATE CrSect_en_elast_m2 !!!")', Rank_of_process
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
        
        DO j = 1, N_en_elast
! was      READ (9, '(4x,f8.2,3x,e9.2)') Energy_en_elast_eV(j), CrSect_en_elast_m2(j)
!!           READ (9, '(4x,f9.3,2x,e10.3)') Energy_en_elast_eV(j), CrSect_en_elast_m2(j)

            read (9,*) Energy_en_elast_eV(j), CrSect_en_elast_m2(j)
!            CrSect_en_elast_m2(j) =  CrSect_en_elast_m2(j) * 1.e-20 !** for Hayashi data

!           print '(4x,f8.2,3x,e9.2)', Energy_en_elast_eV(j), CrSect_en_elast_m2(j)
        END DO

        CLOSE (9, STATUS = 'KEEP')
     
     ELSE

        PRINT '(/2x,"Process ",i3," : Error! Cannot find file ssc_crsect_en_elast.dat")', Rank_of_process      
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
        
     END IF
  
  END IF                                        

! read cross sections of e-n excitation collisions from data file, if this kind of collisions is activated   
  IF (Colflag_kind_spec(2, 1).EQ.1) THEN

     INQUIRE (FILE = 'ssc_crsect_en_excit.dat', EXIST = exists)
     
     IF (exists) THEN
     
        OPEN (9, FILE = 'ssc_crsect_en_excit.dat')
!!        OPEN (19, FILE = 'ssc_crsect_en_excitNIST.dat',status='replace')

        PRINT '(2x,"Process ",i3," : e-n excitation model 1 collisions cross-sections data file is found. Reading the data file...")',&
                                                                                                               & Rank_of_process
        READ (9, '(2x,i4)') N_en_excit
!!        write(19,'(2x,i4)') N_en_excit

        ALLOCATE(Energy_en_excit_eV(1:N_en_excit), STAT=ALLOC_ERR)
        IF(ALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in ALLOCATE Energy_en_excit_eV !!!")', Rank_of_process
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF

        ALLOCATE(CrSect_en_excit_m2(1:N_en_excit), STAT=ALLOC_ERR)
        IF(ALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in ALLOCATE CrSect_en_excit_m2 !!!")', Rank_of_process
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
        
        DO j = 1, N_en_excit
           read(9,*) Energy_en_excit_eV(j), CrSect_en_excit_m2(j)
        END DO
        
        CLOSE (9, STATUS = 'KEEP')
!!        CLOSE (19, STATUS = 'KEEP')       
        DO j = 1, N_en_excit                                 !
           IF (CrSect_en_excit_m2(j).GT.0.0) THEN            !
              Thresh_en_excit_eV = Energy_en_excit_eV(j)     ! finding the energy threshold value
              EXIT                                           !
           END IF                                            !
        END DO                                               !
     
     ELSE

        PRINT '(/2x,"Process ",i3," : Error! Cannot find file ssc_crsect_en_excit.dat !!!")', Rank_of_process      
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
        
     END IF
  
  END IF

! read cross sections of e-n ionization collisions from data file, if this kind of collisions is activated   
  IF (Colflag_kind_spec(3, 1).EQ.1) THEN

     INQUIRE (FILE = 'ssc_crsect_en_ioniz.dat', EXIST = exists)
     IF (exists) THEN
     
        OPEN (9, FILE = 'ssc_crsect_en_ioniz.dat')

        PRINT '(2x,"Process ",i3," : e-n ionization collisions cross-sections data file is found. Reading the data file...")',&
                                                                                                              &  Rank_of_process
        READ (9, '(2x,i4)') N_en_ioniz

        ALLOCATE(Energy_en_ioniz_eV(1:N_en_ioniz), STAT=ALLOC_ERR)
        IF(ALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in ALLOCATE Energy_en_ioniz_eV !!!")', Rank_of_process
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF

        ALLOCATE(CrSect_en_ioniz_m2(1:N_en_ioniz), STAT=ALLOC_ERR)
        IF(ALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in ALLOCATE CrSect_en_ioniz_m2 !!!")', Rank_of_process
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
        
        DO j = 1, N_en_ioniz
! was      READ (9, '(4x,f7.1,4x,e10.3)') Energy_en_ioniz_eV(j), CrSect_en_ioniz_m2(j)
!!           READ (9, '(4x,f8.2,3x,e10.3)') Energy_en_ioniz_eV(j), CrSect_en_ioniz_m2(j)
          read (9,*) Energy_en_ioniz_eV(j), CrSect_en_ioniz_m2(j)
        END DO
        
        CLOSE (9, STATUS = 'KEEP')
        
        DO j = 1, N_en_ioniz                                   !
           IF (CrSect_en_ioniz_m2(j).GT.0.0) THEN              ! 
              Thresh_en_ioniz_eV = Energy_en_ioniz_eV(j)       ! find the energy threshold value
              EXIT                                             !
           END IF                                              !
        END DO                                                 ! 
     
     ELSE

        PRINT '(/2x,"Process ",i3," : Error! Cannot find file ssc_crsect_en_ioniz.dat !!!")', Rank_of_process      
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
        
     END IF
  
  END IF

  ! read cross sections for an additional excitation collision from data file, if this kind of collisions is activated   
  IF (Colflag_kind_spec(5, 1).EQ.1) THEN

     INQUIRE (FILE = 'ssc_crsect_en_excit_2.dat', EXIST = exists)
     IF (exists) THEN
   
        OPEN (9, FILE = 'ssc_crsect_en_excit_2.dat')

        PRINT '(2x,"Process ",i3," : e-n excitation model 2 collisions cross-sections data file is found. Reading the data file...")',&
                                                                                                            &  Rank_of_process
        READ (9, '(2x,i4)') N_en_excit_2

        ALLOCATE(Energy_en_excit_2_eV(1:N_en_excit_2), STAT=ALLOC_ERR)
        IF(ALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in ALLOCATE Energy_en_excit_2_eV !!!")', Rank_of_process
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF

        ALLOCATE(CrSect_en_excit_2_m2(1:N_en_excit_2), STAT=ALLOC_ERR)
        IF(ALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in ALLOCATE CrSect_en_excit_2_m2 !!!")', Rank_of_process
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
      
        DO j = 1, N_en_excit_2
          read (9,*) Energy_en_excit_2_eV(j), CrSect_en_excit_2_m2(j)
        END DO
      
        CLOSE (9, STATUS = 'KEEP')
        DO j = 1, N_en_excit_2                                 !
           IF (CrSect_en_excit_2_m2(j).GT.0.0) THEN            !
              Thresh_en_excit_2_eV = Energy_en_excit_2_eV(j)     ! finding the energy threshold value
              EXIT                                           !
           END IF                                            !
        END DO                                               !

     ELSE

        PRINT '(/2x,"Process ",i3," : Error! Cannot find file ssc_crsect_en_excit_2.dat !!!")', Rank_of_process      
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
      
     END IF

  END IF

END SUBROUTINE CONFIG_READ_CRSECT_ARRAYS 

!=========================================================================================================
SUBROUTINE REMOVE_CRSECT_ARRAYS

  USE ParallelOperationValues
  USE MCCollisions
  IMPLICIT NONE
  INTEGER DEALLOC_ERR

  IF (ALLOCATED(Energy_en_elast_eV)) THEN
     DEALLOCATE(Energy_en_elast_eV, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE Energy_en_elast_eV !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  IF (ALLOCATED(Energy_en_excit_eV)) THEN
     DEALLOCATE(Energy_en_excit_eV, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE Energy_en_excit_eV !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  IF (ALLOCATED(Energy_en_excit_2_eV)) THEN
     DEALLOCATE(Energy_en_excit_2_eV, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE Energy_en_excit_2_eV !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  IF (ALLOCATED(Energy_en_ioniz_eV)) THEN
     DEALLOCATE(Energy_en_ioniz_eV, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE Energy_en_ioniz_eV !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  IF (ALLOCATED(CrSect_en_elast_m2)) THEN
     DEALLOCATE(CrSect_en_elast_m2, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE CrSect_en_elast_m2 !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  IF (ALLOCATED(CrSect_en_excit_m2)) THEN
     DEALLOCATE(CrSect_en_excit_m2, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE CrSect_en_excit_m2 !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  IF (ALLOCATED(CrSect_en_excit_2_m2)) THEN
     DEALLOCATE(CrSect_en_excit_2_m2, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE CrSect_en_excit_2_m2 !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  IF (ALLOCATED(CrSect_en_ioniz_m2)) THEN
     DEALLOCATE(CrSect_en_ioniz_m2, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE CrSect_en_ioniz_m2 !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

END SUBROUTINE REMOVE_CRSECT_ARRAYS

!=========================================================================================================
! Note, this subroutine is called in the INITIATE_MC_COLLISIONS, before calling CONFIGURE_COLLISION_ARRAYS 
!
SUBROUTINE CONFIGURE_COLLSPEC_ARRAYS 

  USE ParallelOperationValues
  USE MCCollisions
  USE CurrentProblemValues, ONLY : N_spec
  IMPLICIT NONE
  INTEGER ALLOC_ERR

!  PRINT '(4x,"Allocating species arrays for electron-neutral and ion-neutral collisions ...")'

  IF (.NOT.ALLOCATED(Ncolkind_spec)) THEN
     ALLOCATE(Ncolkind_spec(1:N_spec), STAT=ALLOC_ERR)
     IF(ALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in ALLOCATE Ncolkind_spec !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  IF (.NOT.ALLOCATED(deltaEnergy_eV_spec)) THEN
     ALLOCATE(deltaEnergy_eV_spec(1:N_spec), STAT=ALLOC_ERR)
     IF(ALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in ALLOCATE deltaEnergy_eV_spec !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

!  PRINT '(2x,"Finished :)")'

END SUBROUTINE CONFIGURE_COLLSPEC_ARRAYS 

!=========================================================================================================
! Note, this subroutine is called in the INITIATE_MC_COLLISIONS, before calling CONFIGURE_COLLISION_ARRAYS 
!
SUBROUTINE SETVALUES_COLLSPEC_ARRAYS 

  USE MCCollisions
  USE CurrentProblemValues, ONLY : N_spec, T_e_eV

  IMPLICIT NONE
!  INTEGER ALLOC_ERR
  INTEGER s          ! species
!  INTEGER i          ! group
  INTEGER k          ! kind of collision

!  PRINT '(4x,"Allocating species arrays for electron-neutral and ion-neutral collisions ...")'

  DO s = 1, N_spec
     deltaEnergy_eV_spec(s) = maxEnergy_eV_spec(s) / (Nenergyval_spec(s) - 1)
  END DO

! determine the number of activated kinds of collisions for different species  
  Ncolkind_spec = 0                                   
  DO s = 1, N_spec                                   
     DO k = 1, maxNcolkind_spec(s)                    
        IF (Colflag_kind_spec(k, s).EQ.1) THEN       
           Ncolkind_spec(s) = Ncolkind_spec(s) + 1   
        END IF                                        
     END DO                                       
  END DO

!  PRINT '(2x,"Finished :)")'

END SUBROUTINE SETVALUES_COLLSPEC_ARRAYS 

!===================================================================================================
! 
SUBROUTINE REMOVE_COLLSPEC_ARRAYS 

  USE ParallelOperationValues
  USE MCCollisions
  USE CurrentProblemValues, ONLY: seed
  IMPLICIT NONE
  INTEGER DEALLOC_ERR

!  PRINT '(4x,"Allocating species arrays for electron-neutral and ion-neutral collisions ...")'

  IF (ALLOCATED(Ncolkind_spec)) THEN
     DEALLOCATE(Ncolkind_spec, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE Ncolkind_spec !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  IF (ALLOCATED(deltaEnergy_eV_spec)) THEN
     DEALLOCATE(deltaEnergy_eV_spec, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE deltaEnergy_eV_spec !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

!  PRINT '(2x,"Finished :)")'

  IF (ALLOCATED(seed)) THEN
     DEALLOCATE(seed, STAT=DEALLOC_ERR)
  END IF

END SUBROUTINE REMOVE_COLLSPEC_ARRAYS 

!=======================================================================================================
! Note, this subroutine is called in the INITIATE_MC_COLLISIONS, after calling SETVALUES_COLLSPEC_ARRAYS  
!
SUBROUTINE CONFIGURE_COLLISION_ARRAYS 

  USE ParallelOperationValues
  USE MCCollisions
  USE CurrentProblemValues, ONLY : N_spec
  IMPLICIT NONE
  INTEGER ALLOC_ERR
  INTEGER s           ! species

!  PRINT '(4x,"Allocating arrays for electron-neutral and ion-neutral collisions ...")'

  IF (.NOT.ALLOCATED(maxColfrac_spec)) THEN
     ALLOCATE(maxColfrac_spec(1:N_spec), STAT=ALLOC_ERR)
     IF(ALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in ALLOCATE maxColfrac_spec !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  IF (.NOT.ALLOCATED(Colkind_of_spec)) THEN
     ALLOCATE(Colkind_of_spec(1:N_spec), STAT=ALLOC_ERR)
     IF(ALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in ALLOCATE Colkind_of_spec !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  DO s = 1, N_spec
     NULLIFY(Colkind_of_spec(s)%activated)
     IF (.NOT.ASSOCIATED(Colkind_of_spec(s)%activated)) THEN
        ALLOCATE(Colkind_of_spec(s)%activated(1:Ncolkind_spec(s)), STAT=ALLOC_ERR)
        IF(ALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in ALLOCATE Colkind_of_spec(s)%activated !!! s = ",i1)', Rank_of_process, s
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
     END IF
  END DO

  IF (.NOT.ALLOCATED(Colprob_of_spec)) THEN
     ALLOCATE(Colprob_of_spec(1:N_spec), STAT=ALLOC_ERR)
     IF(ALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in ALLOCATE Colprob_of_spec !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  DO s = 1, N_spec
     NULLIFY(Colprob_of_spec(s)%kind_energy)
     IF (.NOT.ASSOCIATED(Colprob_of_spec(s)%kind_energy)) THEN
        ALLOCATE(Colprob_of_spec(s)%kind_energy(1:Ncolkind_spec(s),1:Nenergyval_spec(s)), STAT=ALLOC_ERR)
        IF(ALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in ALLOCATE Colprob_of_spec(s)%kind_energy !!! s = ",i1)', Rank_of_process, s
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
     END IF
  END DO

  IF (.NOT.ALLOCATED(Numbers_of_particles)) THEN
     ALLOCATE(Numbers_of_particles(0:(N_of_processes-1)), STAT=ALLOC_ERR)      !^%
     IF (ALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in ALLOCATE Numbers_of_particles !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  IF (.NOT.ALLOCATED(Numbers_of_collisions)) THEN
     ALLOCATE(Numbers_of_collisions(0:(N_of_processes-1)), STAT=ALLOC_ERR)     !^%
     IF (ALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in ALLOCATE Numbers_of_collisions !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF
  

!  PRINT '(2x,"Finished :)")'

END SUBROUTINE CONFIGURE_COLLISION_ARRAYS 
 
!===================================================================================================
! Note, this subroutine is called in the INITIATE_PARAMETERS  
SUBROUTINE SETVALUES_COLLISION_ARRAYS 

  USE ParallelOperationValues
  USE MCCollisions
  USE CurrentProblemValues, ONLY : N_spec, delta_t_s
  IMPLICIT NONE

  REAL(8) FREQUENCY_OF_COLLISION, sum

  INTEGER s   ! species
  INTEGER k   ! index for activated collisions
  INTEGER j   ! index for all collisions, index for energy 
  REAL energy_eV, max_f, aa, bb !, delta_t

!  PRINT '(4x,"Setting the values of arrays for electron-neutral and ion-neutral collisions ...")'

  DO s = 1, N_spec    

! calculate the array of indexes of active collisions Colkind_of_spec(1:N_spec)%activated(1:Ncolkind_spec(s))
     k = 0
     DO j = 1, maxNcolkind_spec(s)              
        IF(Colflag_kind_spec(j, s).EQ.1) THEN
           k = k + 1              
           IF (k.GT.Ncolkind_spec(s)) THEN
              PRINT '(2x,"Process ",i3," : Error in determining Colkind_of_spec(",i2,")%activated")', Rank_of_process, s
              PRINT '(2x,"Current array index ",i2," exceeds the upper limit ",i2)', k, Ncolkind_spec(s)
              PRINT '(2x,"Program will be terminated now :(")'
              STOP
           END IF
           Colkind_of_spec(s)%activated(k) = j     
        END IF                
     END DO  
! if the species has no activated collisions, take the next species
     IF (Ncolkind_spec(s).EQ.0) CYCLE

! calculate the array of probabilities Colprob_of_spec(1:N_spec)%kind_energy(1:Ncolkind_spec(s),1:Nenergyval_spec(s))
     max_f = 0.0_8
     DO j = 1, Nenergyval_spec(s)
        energy_eV = (j-1) * deltaEnergy_eV_spec(s)
        DO k = 1, Ncolkind_spec(s)
      
           Colprob_of_spec(s)%kind_energy(k, j) = FREQUENCY_OF_COLLISION(energy_eV, Colkind_of_spec(s)%activated(k), s)
           IF (k.GT.1) THEN 
              Colprob_of_spec(s)%kind_energy(k, j) = Colprob_of_spec(s)%kind_energy(k, j) + & 
                                                   & Colprob_of_spec(s)%kind_energy(k-1, j)
           END IF

        END DO         !cumulative frequency cycle
        IF (Colprob_of_spec(s)%kind_energy(Ncolkind_spec(s), j).GT.max_f) THEN
          max_f = Colprob_of_spec(s)%kind_energy(Ncolkind_spec(s), j)
        END IF
     END DO

     maxColfrac_spec(s) = max_f
! renormalize the probability (must be not greater than 1), bypassed if every particle collides
!     DO j = 1, Nenergyval_spec(s)
!        DO k = 1, Ncolkind_spec(s)
!           aa = Colprob_of_spec(s)%kind_energy(k, j)
!           Colprob_of_spec(s)%kind_energy(k, j) = aa /  maxColfrac_spec(s)
!        END DO
!!!!!!!!print '(2x,i4,2x,6(1x,e12.5))', j, Colprob_of_spec(s)%kind_energy(1:Ncolkind_spec(s), j)
!     END DO

! calculating the array of fractions of particles, colliding each time step (including the NULL collisions), 
! for different groups and species maxColfrac_spec(1:N_spec)
  
 END DO !species cycle

!  PRINT '(2x,"Finished :)")'

  IF (Rank_of_process.EQ.0) THEN                                              ! server node saves the collision frequencies

     IF (Colflag_kind_spec(1, 1).EQ.1) THEN                                   ! electron-neutral, elastic
        OPEN (99, FILE = '_en_elast_coll_freqs.dat')
        DO j = 1, Nenergyval_spec(1)
           energy_eV = (j-1) * deltaEnergy_eV_spec(1)
           WRITE (99, '(2(2x,e12.5))') energy_eV, FREQUENCY_OF_COLLISION(energy_eV, 1, 1) / delta_t_s
        END DO
        CLOSE (99, STATUS = 'KEEP')
     END IF

     IF (Colflag_kind_spec(2, 1).EQ.1) THEN                                   ! electron-neutral, excitation model 1
        OPEN (99, FILE = '_en_excit_coll_freqs.dat')
        DO j = 1, Nenergyval_spec(1)
           energy_eV = (j-1) * deltaEnergy_eV_spec(1)
           WRITE (99, '(2(2x,e12.5))') energy_eV, FREQUENCY_OF_COLLISION(energy_eV, 2, 1) / delta_t_s
        END DO
        CLOSE (99, STATUS = 'KEEP')
     END IF

     IF (Colflag_kind_spec(3, 1).EQ.1) THEN                                   ! electron-neutral, ionization
        OPEN (99, FILE = '_en_ioniz_coll_freqs.dat')
        DO j = 1, Nenergyval_spec(1)
           energy_eV = (j-1) * deltaEnergy_eV_spec(1)
           WRITE (99, '(2(2x,e12.5))') energy_eV, FREQUENCY_OF_COLLISION(energy_eV, 3, 1) / delta_t_s
!print '(2(2x,e12.5))', energy_eV, FREQUENCY_OF_COLLISION(energy_eV, 3, 1) / delta_t_s 
        END DO
        CLOSE (99, STATUS = 'KEEP')
     END IF

     IF (Colflag_kind_spec(4, 1).EQ.1) THEN                                   ! electron, turbulent
        OPEN (99, FILE = '_e_turb_coll_freqs.dat')
        DO j = 1, Nenergyval_spec(1)
           energy_eV = (j-1) * deltaEnergy_eV_spec(1)
           WRITE (99, '(2(2x,e12.5))') energy_eV, FREQUENCY_OF_COLLISION(energy_eV, 4, 1) / delta_t_s
        END DO
        CLOSE (99, STATUS = 'KEEP')
      END IF

     IF (Colflag_kind_spec(4, 1).EQ.1) THEN                                   ! electron, excitation model 2
        OPEN (99, FILE = '_e_excit_2_coll_freqs.dat')
        DO j = 1, Nenergyval_spec(1)
           energy_eV = (j-1) * deltaEnergy_eV_spec(1)
           WRITE (99, '(2(2x,e12.5))') energy_eV, FREQUENCY_OF_COLLISION(energy_eV, 5, 1) / delta_t_s
        END DO
        CLOSE (99, STATUS = 'KEEP')
      END IF

     IF (Colflag_kind_spec(2, 2).EQ.1) THEN                                  ! ion, charge exchange 
       OPEN (99, FILE = '_in_chargex_coll_freqs.dat')
       DO j = 1, Nenergyval_spec(2)
          energy_eV = (j-1) * deltaEnergy_eV_spec(2)
          WRITE (99, '(2(2x,e12.5))') energy_eV, FREQUENCY_OF_COLLISION(energy_eV, 2, 2) / delta_t_s
       END DO
       CLOSE (99, STATUS = 'KEEP')
     END IF

     OPEN (99, FILE = '_en_cumul_probab.dat')
     DO j = 1, Nenergyval_spec(1)
       energy_eV = (j-1) * deltaEnergy_eV_spec(1)
       WRITE (99, '(4(2x,e12.5))') energy_eV, (Colprob_of_spec(1)%kind_energy(k,j), k=1,Ncolkind_spec(1))
       END DO
     CLOSE (99, STATUS = 'KEEP')

     OPEN (99, FILE = '_in_cumul_probab.dat')
     DO j = 1, Nenergyval_spec(2)
       energy_eV = (j-1) * deltaEnergy_eV_spec(2)
       WRITE (99, '(4(2x,e12.5))') energy_eV, (Colprob_of_spec(2)%kind_energy(k,j), k=1,Ncolkind_spec(2))
       END DO
     CLOSE (99, STATUS = 'KEEP') 
     

  END IF

END SUBROUTINE SETVALUES_COLLISION_ARRAYS 

!===================================================================================================
SUBROUTINE REMOVE_COLLISION_ARRAYS

  USE ParallelOperationValues
  USE MCCollisions
  USE CurrentProblemValues, ONLY : N_spec
  IMPLICIT NONE
  INTEGER DEALLOC_ERR
  INTEGER s    ! species

!  PRINT '("Deleting arrays for electron-neutral and ion-neutral collisions ...")'

  IF (ALLOCATED(maxColfrac_spec)) THEN
     DEALLOCATE(maxColfrac_spec, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE maxColfrac_spec !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  DO s = 1, N_spec
     IF (ASSOCIATED(Colkind_of_spec(s)%activated)) THEN
        DEALLOCATE(Colkind_of_spec(s)%activated, STAT=DEALLOC_ERR)
        IF(DEALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE Colkind_of_spec(s)%activated !!! s = ",i1)', Rank_of_process, s
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
     END IF
  END DO

  IF (ALLOCATED(Colkind_of_spec)) THEN
     DEALLOCATE(Colkind_of_spec, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE Colkind_of_spec !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  DO s = 1, N_spec
     IF (ASSOCIATED(Colprob_of_spec(s)%kind_energy)) THEN
        DEALLOCATE(Colprob_of_spec(s)%kind_energy, STAT=DEALLOC_ERR)
        IF(DEALLOC_ERR.NE.0)THEN
           PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE Colprob_of_spec(s)%kind_energy !!! s = ",i1)', Rank_of_process, s
           PRINT  '(2x,"Program will be terminated now :(")'
           STOP
        END IF
     END IF
  END DO

  IF (ALLOCATED(Colprob_of_spec)) THEN
     DEALLOCATE(Colprob_of_spec, STAT=DEALLOC_ERR)
     IF(DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE Colprob_of_spec !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF


  IF (ALLOCATED(Numbers_of_particles)) THEN
     DEALLOCATE(Numbers_of_particles, STAT=DEALLOC_ERR)
     IF (DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE Numbers_of_particles !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  IF (ALLOCATED(Numbers_of_collisions)) THEN
     DEALLOCATE(Numbers_of_collisions, STAT=DEALLOC_ERR)
     IF (DEALLOC_ERR.NE.0)THEN
        PRINT '(/2x,"Process ",i3," : Error in DEALLOCATE Numbers_of_collisions !!!")', Rank_of_process
        PRINT  '(2x,"Program will be terminated now :(")'
        STOP
     END IF
  END IF

  NULLIFY(Collided_particle)   

!  PRINT '(2x,"Finished :)")'

END SUBROUTINE REMOVE_COLLISION_ARRAYS

!==========================================================
! 
SUBROUTINE PROCESS_COLLISIONS_WITH_ATOMS

  use mpi
  USE ParallelOperationValues
  USE MCCollisions
  USE CurrentProblemValues
  USE mt19937
  IMPLICIT NONE

!  INCLUDE 'mpif.h'

  INTERFACE
     RECURSIVE SUBROUTINE Node_Killer(node)
       USE MCCollisions
       TYPE (binary_tree), POINTER :: node
     END SUBROUTINE Node_Killer
  END INTERFACE

  REAL RAN, RANF

  INTEGER s             ! species 
  INTEGER random_j      ! number of selected particle
  REAL    real_random_j ! real value used to avoid errors if the number becomes more than the largest 32-bit integer 
  INTEGER j             ! counter of attempts to collide

  REAL    R_collided ! number of collided particles during one timestep, can be less than 1
  INTEGER I_collided ! integer part of R_collided
  REAL    F_collided ! fractional part of R_collided

! one should definitely take I_collided particles as the collided, and 
! one particle with the probability F_collided

  INTEGER start_attempt ! start index in the cycle, which takes the random numbers of particles
                        ! it is 1 if F_collided = 0 and it is 0 if F_collided > 0

  REAL    random_r, R   ! random probability
  INTEGER k          ! counter of kinds of collisions 
 
  REAL(8) energy_eV     ! energy of selected particle [eV]

  INTEGER indx_energy   ! 2-nd index of the array Colprob_of_spec(s)%kind_energy(k, indx_energy)  
  INTEGER indx_coll     ! index of selected collision type

  REAL(8) Vx_n, Vy_n, Vz_n  ! velociity components of the colliding neutral atom

  LOGICAL search_again
  LOGICAL already_stored
  INTEGER ALLOC_ERR
  LOGICAL Find_in_stored_list

  INTEGER n, isum                           ! process number
  INTEGER ierr

  INTEGER N_of_collisions        ! number of collisions for a client process

! function
  REAL(8) Shape_of_Natom  ! defines shape of the neutral atom density profile, 0<=Shape_of_Natom<=1
  
  call random_seed(put=seed)

  DO s = 1, N_spec

     IF (Ncolkind_spec(s).EQ.0) CYCLE                  ! skip if no collision for the species

     N_of_collisions = 0

     IF (Rank_of_process.EQ.0) THEN                                                              ! server >>>

        Numbers_of_particles = 0
        N_part = 0
! Receive data from clients
        CALL MPI_GATHER(N_part(s), 1, MPI_INTEGER, Numbers_of_particles, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
        
        DO n = 2, N_of_processes - 1        
           Numbers_of_particles(n) = Numbers_of_particles(n) + Numbers_of_particles(n - 1) 
        END DO
        N_part(s) = Numbers_of_particles(N_of_processes - 1)      ! ??? was (n-1)

        R_collided = maxColfrac_spec(s) * N_part(s)      
        I_collided = int(R_collided)
        F_collided = R_collided - I_collided

        IF (F_collided.GT.0.000001) THEN
           start_attempt = 0                 ! start with statistical collisions
        ELSE
           start_attempt = 1                 ! skip the statistical collisions if the fractional part is too small (to avoid overflow)
        END IF

!============= Randomly take some particle "I_collided(+1)" times, define number of collisions for each process
        Numbers_of_collisions = 0
        DO j = start_attempt, I_collided

           IF (j.EQ.0) THEN                                                     ! statistical part
              real_random_j = RAN(I_random_seed) * N_part(s) / F_collided       ! the value of the right-hand part can exceed the maximal 32-bit integer number
              IF ((real_random_j.GT.N_part(s)).OR.(real_random_j.LT.0.0)) CYCLE ! that is why we use the real value and check the negative value here too
              random_j = real_random_j                                          ! initialize the integer value
              IF (random_j.GT.N_part(s)) CYCLE                                  ! double check (maybe that's silly...)
           ELSE                                                                 ! mandatory part
              random_j = RAN(I_random_seed) * N_part(s)                         ! note, we assume that N_part(s) is a trusted 32-bit integer value
           END IF                                                           
           IF (random_j.LT.1) random_j = 1                                      ! to be sure that we are within the array boundaries
           
           DO n = 1, N_of_processes - 1
              IF (random_j.LE.Numbers_of_particles(n)) THEN 
                 Numbers_of_collisions(n) = Numbers_of_collisions(n) + 1
                 EXIT
              END IF
           END DO

        END DO

! TRANSMIT DATA TO THE CLIENTS
        CALL MPI_SCATTER(Numbers_of_collisions, 1, MPI_INTEGER, N_of_collisions, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

     ELSE                                                   ! client >>>

! send data to the server
        CALL MPI_GATHER(N_part(s), 1, MPI_INTEGER, Numbers_of_particles, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

! Allocate binary tree to store the numbers of particles which have collided already
        NULLIFY(Collided_particle)
        IF (.NOT.ASSOCIATED(Collided_particle)) THEN
           ALLOCATE(Collided_particle, STAT=ALLOC_ERR)
           IF (ALLOC_ERR.NE.0) THEN
              PRINT '(/2x,"Process ",i3," : Error in ALLOCATE Collided_particle !!!")', Rank_of_process
              PRINT  '(2x,"Program will be terminated now :(")'
              STOP
           END IF
        END IF
        NULLIFY(Collided_particle%Larger)
        NULLIFY(Collided_particle%Smaller)

! receive data from the server
        CALL MPI_SCATTER(Numbers_of_collisions, 1, MPI_INTEGER, N_of_collisions, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

!============= Randomly take some (always different) particle "I_collided(+1)" times
!        DO j = 1, N_of_collisions
          DO j = 1, N_part(s) ! go over all particles instead
            
            random_j = j !not random here, all particles collide

!           search_again = .TRUE.                              ! if I_collided >= 1 then
!           DO WHILE (search_again)                            ! search will be repeated until a number will be successfully obtained
!              random_j = int(RAN(I_random_seed) * N_part(s))
!              IF (random_j.LT.1) random_j = 1                !                   
!              IF (random_j.GT.N_part(s)) CYCLE                ! within the array boundaries
!              already_stored = Find_in_stored_list(random_j)   ! 
!              IF (already_stored) CYCLE                                  ! skip particle if it has collided already  
!              search_again = .FALSE.
!           END DO
           
!------------- Determine the kind of collision for the selected particle
!!           random_r = RAN(I_random_seed)
!!           random_r = RANF()
!!           call random_number(random_r)
           random_r = grnd()

           SELECT CASE (s)
           CASE (1)

! #### NEW: first check whether collision is allowed by the neutral density profile #####
!!              call random_number(R)
              R = grnd()
              IF (DBLE(R).LE.Shape_of_Natom(X_of_spec(s)%part(random_j))) THEN      
! if we are here, we can proceed with processing the colliding particle as usual           
! for electrons the kind of collision is determined by the energy of electron
                 energy_eV = energy_factor_eV_s(s) * (VX_of_spec(s)%part(random_j)**2 + &
                                                   &  VY_of_spec(s)%part(random_j)**2 + & 
                                                   &  VZ_of_spec(s)%part(random_j)**2)

                 indx_energy = int(energy_eV / deltaEnergy_eV_spec(s)) + 1
                 IF (indx_energy.GT.Nenergyval_spec(s)) THEN
                    indx_energy = Nenergyval_spec(s)
                    !                    PRINT *, '(2x,"particle of species ",i2,"exceeds the energy limit for collisions")', s
                 END IF
              
                 DO k = Ncolkind_spec(s), 1, -1 
                    IF (random_r.GT.Colprob_of_spec(s)%kind_energy(k, indx_energy)) EXIT 
                 END DO
              
                 indx_coll = k + 1        
                 CALL COLLIDE_ELECTRON(indx_coll, random_j, energy_eV)
              END IF

           CASE (2)   
! for ions the kind of collision is determined by the energy of relative motion ion-neutral
! that's why we must take the neutral's velocity from Maxwell distribution:
              CALL GetMaxwellVelocity(Vx_n) 
              CALL GetMaxwellVelocity(Vy_n) 
              CALL GetMaxwellVelocity(Vz_n) 

!!              CALL GetMaxwellVector(Vx_n, Vy_n, Vz_n)

! Use the factor above to obtain the dim-less velocity (V/V_te) of the produced neutral atom
              Vx_n = Vx_n * alpha_Vscl
              Vy_n = Vy_n * alpha_Vscl
              Vz_n = Vz_n * alpha_Vscl
! Take the energy of ion at the frame where the neutral atom is at rest
! Ngas * sigma(|v-u|) * |v-u| will be looked up in the collision frequency array
              energy_eV = energy_factor_eV_s(s) * ( (VX_of_spec(s)%part(random_j) - Vx_n)**2 + &
                                                &   (VY_of_spec(s)%part(random_j) - Vy_n)**2 + & 
                                                &   (VZ_of_spec(s)%part(random_j) - Vz_n)**2 )

              indx_energy = int(energy_eV / deltaEnergy_eV_spec(s)) + 1
              IF (indx_energy.GT.Nenergyval_spec(s)) THEN
                 indx_energy = Nenergyval_spec(s)
!                    PRINT *, '(2x,"particle of species ",i2,"exceeds the energy limit for collisions")', s
              END IF
! *** the type of collision is chosen here:
              DO k = Ncolkind_spec(s), 1, -1 
                 IF (random_r.GT.Colprob_of_spec(s)%kind_energy(k, indx_energy)) EXIT 
              END DO
           
              indx_coll = k + 1
              CALL COLLIDE_ION(indx_coll, random_j, Vx_n, Vy_n, Vz_n) 
!          CASE (3)
           END SELECT    !------------- End of determination of the kind of collision

        END DO           !============= End of random particle taking

        CALL Node_Killer(Collided_particle) !*** null-collision bypassed

     END IF !choice between species 

  END DO                 !------------- End of the external cycle over species

!print *, 'exit  PROCESS_COLLISIONS_WITH_ATOMS' ------------------------------------------------------------------------------------
  call random_seed(get=seed)
  return
END SUBROUTINE PROCESS_COLLISIONS_WITH_ATOMS

!-----------------------------------------------------------------
! function's value equals .TRUE. if the particle is already stored,
! otherwise function's value is .FALSE. (i.e. particle does not collide yet)
LOGICAL FUNCTION Find_in_stored_list(number)

  USE MCCollisions
  IMPLICIT NONE

  INTEGER number
  TYPE (binary_tree), POINTER :: current

  Find_in_stored_list = .FALSE.

  current => Collided_particle

  DO 

     IF (number.GT.current%number) THEN
        IF (ASSOCIATED(current%Larger)) THEN
           current => current%Larger               ! 
           CYCLE                                   ! go to the next node, with larger "number"
        ELSE
           EXIT
        END IF
     END IF

     IF (number.LT.current%number) THEN
        IF (ASSOCIATED(current%Smaller)) THEN
           current => current%Smaller              ! 
           CYCLE                                   ! go to the next node, with smaller "number"
        ELSE
           EXIT
        END IF
     END IF

     Find_in_stored_list = .TRUE.                  ! number.EQ.current%number
     EXIT                                          ! if we are here, then we found the match
     
  END DO

END FUNCTION Find_in_stored_list

!-----------------------------------------------------------------
! subroutine adds number to the binary tree
! we assume that there are no nodes in the tree with the same value yet
SUBROUTINE Add_to_stored_list(number)

  USE ParallelOperationValues
  USE MCCollisions

  IMPLICIT NONE
  INTEGER number
  TYPE (binary_tree), POINTER :: current
  INTEGER ALLOC_ERR

  current => Collided_particle                  ! start from the head node of the binary tree

  DO                                            ! go through the allocated nodes to the end of the branch

     IF (number.GT.current%number) THEN         
        IF (ASSOCIATED(current%Larger)) THEN        
           current => current%Larger               ! 
           CYCLE                                   ! go to the next node, with larger "number"
        ELSE
           ALLOCATE(current%Larger, STAT=ALLOC_ERR)
           IF (ALLOC_ERR.NE.0) THEN
              PRINT '(/2x,"Process ",i3," : Add_to_stored_list: Error in ALLOCATE current%Larger !!!")', Rank_of_process
              PRINT  '(2x,"Program will be terminated now :(")'
              STOP
           END IF
           current => current%Larger
           EXIT
        END IF
     END IF

     IF (number.LT.current%number) THEN
        IF (ASSOCIATED(current%Smaller)) THEN        
           current => current%Smaller              ! 
           CYCLE                                   ! go to the next node, with smaller "number"
        ELSE
           ALLOCATE(current%Smaller, STAT=ALLOC_ERR)
           IF (ALLOC_ERR.NE.0) THEN
              PRINT '(/2x,"Process ",i3," : Add_to_stored_list: Error in ALLOCATE current%Smaller !!!")', Rank_of_process
              PRINT  '(2x,"Program will be terminated now :(")'
              STOP
           END IF
           current => current%Smaller
           EXIT
        END IF
     END IF
     
  END DO

  current%number = number                       ! store the number
  NULLIFY(current%Larger)
  NULLIFY(current%Smaller)

END SUBROUTINE Add_to_stored_list

!---------------------------------------------------
! this subroutine kills the nodes of the binary tree
RECURSIVE SUBROUTINE Node_Killer(node)

  USE ParallelOperationValues
  USE MCCollisions
  IMPLICIT NONE
  
  TYPE (binary_tree), POINTER :: node
  INTEGER DEALLOC_ERR

  IF (ASSOCIATED(node%Larger))  CALL Node_Killer(node%Larger)
  IF (ASSOCIATED(node%Smaller)) CALL Node_Killer(node%Smaller)

  DEALLOCATE(node, STAT=DEALLOC_ERR)
  IF (DEALLOC_ERR.NE.0) THEN
     PRINT '(/2x,"Process ",i3," : Error in Node_Killer !!!")', Rank_of_process
     PRINT  '(2x,"Program will be terminated now :(")'
     STOP
  END IF

  RETURN

END SUBROUTINE Node_Killer


!*******************************************************************************************
! This module is used by two subroutines, calculating the arbitrary velocity
! according to the isotropic maxwell distribution
MODULE MaxwellVelocity
  REAL(8) v(0:9000)          !(0:300)
  REAL(8) v_inj(0:4500)      !(0:150)
  REAL(8) :: U_max = 3.0_8
  INTEGER :: R_max     = 9000 !300
  INTEGER :: R_max_inj = 4500 !150
END MODULE MaxwellVelocity

!-------------------------------------------------------------------------------------------
! Prepares the tabulated values of integral of the maxwell distribution function
!  
SUBROUTINE PrepareMaxwellDistribIntegral

  USE ParallelOperationValues
  USE MaxwellVelocity
  USE CurrentProblemValues, ONLY : N_box_vel
  IMPLICIT NONE

  INTEGER i
  INTEGER N_pnts
  INTEGER count
  REAL(8) V_min, V_max
  REAL(8) F(0:180003)  !(0:30003)            ! to be sure that we overcome V_max
  REAL(8) temp
  REAL(8) dV

  LOGICAL check1, check2

  check1 = .FALSE.
  check2 = .FALSE.

! ------- for symmetrical maxwellian
  N_pnts = 180000  !30000
  V_min  = -U_max
  V_max  =  U_max
  dV = (V_max - V_min) / N_pnts

  F = 0.0_8
  DO i = 1, N_pnts + 3
     F(i) = F(i-1) + EXP( - (V_min + (REAL(i)-0.5_8) * dV)**2 )
  END DO

  temp = F(N_pnts)
  F = F * R_max / temp   ! normalize integral such that F(N_pnts) = R_max

  v(0) = V_min
  count = 0
  open(unit=99, file = 'SampledMaxwell.dat', status = 'replace')
  DO i = 1, N_pnts + 3
     IF ((INT(F(i))-count).EQ.1) THEN
        count = count + 1
        v(count) = V_min + i * dV
        write (99,*) count, v(count)
        IF (count.EQ.R_max) THEN
           check1 = .TRUE.
           EXIT
        END IF
     END IF
  END DO
  close(unit=99)
 
  v = v * N_box_vel

!--------- for asymmetrical maxwellian * v (used for injection, v > 0)

  N_pnts = 90000   !15000
  V_min  = 0.0_8
  V_max  = U_max
  dV = (V_max - V_min) / N_pnts

  F = 0.0_8
  DO i = 1, N_pnts + 3
     temp = V_min + (REAL(i)-0.5_8) * dV
     F(i) = F(i-1) + EXP( - temp**2 ) * temp
  END DO

  temp = F(N_pnts)
  F(1:(N_pnts+3)) = F(1:(N_pnts+3)) * R_max_inj / temp   ! normalize integral such that F(N_pnts) = R_max_inj

  v_inj(0) = V_min
  count = 0
  DO i = 1, N_pnts + 3
     IF ((INT(F(i))-count).EQ.1) THEN
        count = count + 1
        v_inj(count) = V_min + i * dV
        IF (count.EQ.R_max_inj) THEN
           check2 = .TRUE.
           EXIT
        END IF
     END IF
  END DO

  v_inj = v_inj * N_box_vel

  IF (check1.AND.check2) THEN
     PRINT '(2x,"Process ",i3," : Integrals for producing maxwell distributions are successfully obtained ...")', &
                                                                                                  & Rank_of_process
  ELSE
     PRINT '(2x,"Process ",i3," : ERROR in PrepareMaxwellDistribIntegral !!!")', Rank_of_process
     PRINT '(2x,"The initialization in PrepareMaxwellDistribIntegral is not performed !!!")'
     PRINT '(2x,"The program will be terminated now :(")'
     STOP
  END IF

END SUBROUTINE PrepareMaxwellDistribIntegral

!-------------------------------------------------------------------------------------------
!  
SUBROUTINE GetMaxwellVelocity(U) 
!! the Box-Muller algorithm
  USE mt19937
  USE CurrentProblemValues, ONLY : seed, N_box_vel 
  IMPLICIT NONE

  logical flag
  REAL(8) U, dR1, dR2, W
!!  call random_seed(put=seed)
  flag = .true.
  do while (flag)
!!    call random_number(R1)
!!    call random_number(R2)
!!    dR1 = 2.0_8 * dble(R1) - 1.0_8
!!    dR2 = 2.0_8 * dble(R2) - 1.0_8
    dR1 = 2.0_8 * grnd() - 1.0_8
    dR2 = 2.0_8 * grnd() - 1.0_8
    W = dR1*dR1 + dR2*dR2
    if (W.le.1.d0 .and. W.gt.0.d0) then
      U = -dR1 / dsqrt(W) * dsqrt( -log(W))
! sampling the Gaussian with sigma^2 = 0.5, in accordance with
! how the thermal velocity is defined in the code
      flag = .false.
    else
    end if
  end do
  U = U * dble(N_box_vel) 
!!  call random_seed(get=seed)
  RETURN
  END SUBROUTINE GetMaxwellVelocity

 SUBROUTINE GetMaxwellVector(VX,VY,VZ)
!!*** needs to be validated

 USE MaxwellVelocity
 USE CurrentProblemValues, ONLY : I_random_seed, N_box_vel
 Implicit none
 Real(8) U, VX, VY, VZ, pi, twopi
 Real RAN, R1, R2, R3, phi
 logical flag

 pi = 3.14159265358979_8
 twopi = 2. * pi
 flag = .true.
 
 do while (flag)
   R1 = RAN(I_random_seed)
   R2 = RAN(I_random_seed)
   R3 = RAN(I_random_seed)
   if (R1.ge.1.e-6 .and. R2.ge.1.e-6) then
     U = -( log(R1) + log (R2) * 0.5 * (1. + cos(pi * R3)) )
     flag = .false.
   else
   end if
 end do

 R1 = 1. - 2. * RAN(I_random_seed)
 phi = twopi * RAN(I_random_seed)
 U = sqrt(U) * dble(N_box_vel)
 VX = U * R1 
 VY = U * sqrt(1.-R1*R1)
 VZ = VY * sin(phi)
 VY = VY * cos(phi)
  
 RETURN
 END SUBROUTINE GetMaxwellVector
!__________________________________________________________________________________________
FUNCTION RANF() RESULT (CR)
!
! Function to generate a uniform random number in [0,1]
! following x(i+1)=a*x(i) mod c with a=7** 5 and
! c=2** 31-1.  Here the seed is a global variable.
!
  USE CurrentProblemValues, ONLY : I_random_seed 
  IMPLICIT NONE
  INTEGER :: H, L, T, A, C, Q, R, SEED
  DATA A/16807/, C/2147483647/, Q/127773/, R/2836/
  REAL :: CR
!
  SEED = I_random_seed
  H = SEED/Q
  L = MOD(SEED, Q)
  T = A*L - R*H
  IF (T .GT. 0) THEN
    SEED = T
  ELSE
    SEED = C + T
  END IF
  CR = SEED/FLOAT(C)
  I_random_seed = SEED
END FUNCTION RANF


!-------------------------------------------------------------------------------------------
!
!  
SUBROUTINE GetInjMaxwellVelocity(U) 

  USE MaxwellVelocity
  USE CurrentProblemValues, ONLY : I_random_seed 
  USE mt19937
  IMPLICIT NONE

  REAL(8) U
  REAL RAN, R
  INTEGER indx
  
!!  R = R_max_inj * RAN(I_random_seed)
  R = R_max_inj * grnd()
  indx = INT(R)

  IF (indx.LT.R_max_inj) THEN
     U = v_inj(indx) + (R - indx) * (v_inj(indx+1) - v_inj(indx))
  ELSE
     U = v_inj(R_max_inj)
  END IF
  RETURN
  
END SUBROUTINE GetInjMaxwellVelocity

!--------------------------------------------------
!
REAL(8) FUNCTION Shape_of_Natom(x)

  USE CurrentProblemValues, ONLY : N_cells, delta_x_m

  IMPLICIT NONE

  REAL(8) x
  
!  Shape_of_Natom = MIN( 400.0_8**(((x-DBLE(N_cells)) * delta_x_m) / 1.0_8), 1.0_8)

!  Shape_of_Natom = MIN( 10.0_8**(((x-DBLE(N_cells)) * delta_x_m) / 0.05_8), 1.0_8)  ! was 400.0_8 instead of 10.0_8

!  Shape_of_Natom = MIN( 4.0_8**(((x-DBLE(N_cells)) * delta_x_m) / 0.1_8), 1.0_8)  ! was 400.0_8 instead of 10.0_8


! here 400 and 1 are because the density decreases 400 times 
! if one goes 1 meter leftward from the right end (N_cells*delta_x_m) 
! also, shape_of_natom should not exceed 1

  Shape_of_Natom = 1.0_8  ! use this in case of uniform neutral density 

END FUNCTION Shape_of_Natom

SUBROUTINE GetMaxwell2D(U, V)
!! the Box-Muller algorithm
  USE mt19937 !** Mercenne twister
  USE CurrentProblemValues, ONLY : N_box_vel

  IMPLICIT NONE
  logical tryagain
  REAL(8) U, V, dR1, dR2, W, aa

  tryagain = .true.
  do while (tryagain)
    dR1 = 2.0_8 * grnd() - 1.0_8
    dR2 = 2.0_8 * grnd() - 1.0_8
    W = dR1*dR1 + dR2*dR2
    if (W.le.1.0_8 .and. W.gt.0.0_8) then
      aa = dsqrt(-dlog(W) / W)
      U = dR1 * aa
      V = dR2 * aa
! sampling the Gaussian with sigma^2 = 0.5, in accordance with
! how the thermal velocity is defined in the code
      tryagain = .false.
    else
    end if
  end do
  U = U * dble(N_box_vel)
  V = V * dble(N_box_vel)
  RETURN
  END SUBROUTINE GetMaxwell2D





