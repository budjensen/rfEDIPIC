FC := mpif90
CC := mpicc
# Get the actual Fortran compiler invoked by mpif90. Need to use double $$ in awk.
BASE_FCOMP := $(shell mpif90 -show | awk '{print $$1}' | awk -F/ '{print $$NF}' | awk -F- '{print $$NF}')

ifeq ($(BASE_FCOMP),ifort)
   $(info "Using Intel compiler" )
   FFLAGS := -O2
   CFLAGS := -O2
else ifeq ($(BASE_FCOMP),gfortran)
   $(info "Using gcc/gfortran compilers" )
   FFLAGS := -O2 -ffree-line-length-none
   CFLAGS := -O2
else
   FFLAGS := -O2
   CFLAGS := -O2
endif

OBJ := parMT19937.o \
       parModules_RF.o \
       parBeamInjection_newflux.o \
       parBeamInPlasma_newflux.o \
       parCheckpoint.o \
       parCollisionProc.o \
       parCurProbVals_dt_NewInt.o \
       parDiagNoGraph_RF.o \
       parFieldCalcProc_RF_diel.o \
       parIonInducedSEEmission_newflux.o \
       parLangevinColl.o \
       parMagneticFieldProfiles_zero.o \
       parMainProgram.o \
       parMCCollisions_BoxMuller_orig.o \
       parPartDynamicProc_Gauss.o \
       parPartExchangeProc.o \
       parPlasmaSourceProc_newflux.o \
       parRandomInit.o \
       parSEEmission_newflux.o \
       parSnapshot_newflux_heat.o

dedipic : $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

%.o : %.f90
	$(FC) $(FFLAGS) -c $< 

%.o : %.F90
	$(FC) $(FFLAGS) -c $< 

%.o : %.c
	$(CC) $(CFLAGS) -c $<

debug: FFLAGS += -g  # Add debugging flags to FFLAGS
debug: CFLAGS += -g  # Add debugging flags to CFLAGS
debug: dedipic       # Build the debug version

clean:
	rm -f *.mod *.o dedipic

