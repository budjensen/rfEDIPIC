FC := mpif90
CC := mpicc
# Get the actual Fortran compiler invoked by mpif90. Need to use double $$ in awk.
BASE_FCOMP := $(shell mpif90 -show | awk '{print $$1}' | awk -F/ '{print $$NF}')

ifeq ($(BASE_FCOMP),ifort)
   $(info "Using Intel compiler" )
   FFLAGS := -O2
   CFLAGS := -O2
else ifeq ($(BASE_FCOMP),gfortran)
   $(info "Using gcc/gfortran compilers" )
   FFLAGS := -O2 -ffree-line-length-none -I${PETSC_DIR}/include
   CFLAGS := -O2
else
   FFLAGS := -O2
   CFLAGS := -O2
endif

OBJ := WELL19937a_new.o \
       parModules.o \
       parRandomNumberInterface.o \
       parBeamInjection.o \
       parBeamInPlasma.o \
       parCheckpoint.o \
       parCollisionProc.o \
       parCurProbVals.o \
       parDiagNoGraph.o \
       parFieldCalcProc.o \
       parIonInducedSEEmission.o \
       parLangevinColl.o \
       parMainProgram.o \
       parMCCollisions.o \
       parPartDynamicProc.o \
       parPartExchangeProc.o \
       parPlasmaSourceProc.o \
       parSEEmission.o \
       parSnapshot.o

dedipic : $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

%.o : %.f90
	$(FC) $(FFLAGS) -c $< 

%.o : %.F90
	$(FC) $(FFLAGS) -c $< 

%.o : %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.mod *.o dedipic

