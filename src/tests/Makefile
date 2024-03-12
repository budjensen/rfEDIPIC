FC := mpif90
SRC_DIR := .
PARENT_DIR := ..

# Get the actual Fortran compiler invoked by mpif90. Need to use double $$ in awk.
BASE_FCOMP := $(shell mpif90 -show | awk '{print $$1}' | awk -F/ '{print $$NF}' | awk -F- '{print $$NF}')

ifeq ($(BASE_FCOMP),ifort)
   $(info "Using Intel compiler" )
   FFLAGS := -O2
else ifeq ($(BASE_FCOMP),gfortran)
   $(info "Using gcc/gfortran compilers" )
   FFLAGS := -O2 -ffree-line-length-none
else
   FFLAGS := -O2
endif

IGNORE_FILES := parMainProgram.f90 parRFDiagnostics.f90 parHeating.f90 parMT19937.f90 parModules.f90

# Specify the order of compilation for first.f90 and second.f90
SRC_FILES := $(PARENT_DIR)/parMT19937.f90 \
             $(PARENT_DIR)/parModules.f90 \
             $(filter-out $(addprefix $(PARENT_DIR)/,$(IGNORE_FILES)), $(wildcard $(PARENT_DIR)/*.f90)) \
             $(wildcard *.f90)

TEMP := $(SRC_FILES:%.f90=%.o)
OBJ_FILES := $(subst $(PARENT_DIR)/,,$(TEMP))


tests: $(OBJ_FILES)
	$(FC) $(FFLAGS) -o $@ $^

%.o: $(PARENT_DIR)/%.f90
	$(FC) $(FFLAGS) -c $<

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

.PHONY: clean
clean:
	rm -f tests *.o *.mod
