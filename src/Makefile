# Paths definition
IDIR = ../include
BDIR = ../bin
ODIR = ./obj
MDIR = ./mod
DIR  = ~/mis_bin
libpath= -lblas -llapack

# Compiler/Linker settings
FC = gfortran -I$(IDIR) -J$(MDIR)  # -Wno-argument-mismatch 
FLFLAGS = -g
FCFLAGS = -g -c -Wall -Wextra -Wconversion -Og -fcheck=bounds -fmax-errors=5
PROGRAM = atomtypes
SRCF90  = equivalent.f90
SRCMOD  = lengths.f90       \
          units.f90         \
          timings.f90       \
          printings.f90     \
          utils.f90         \
          datatypes.f90     \
          sorting.f90       \
          proc_inp.f90      \
          g16_files.f90     \
          gromacs_files.f90 \
          graphtools.f90    \
          isomorphism.f90   \
          genforcefield.f90 \
          wrap_lapack.f90

# Find all source files, create a list of corresponding object files
OBJ     = $(patsubst %.f90,$(ODIR)/%.o,$(SRCF90))
OBJ    := $(OBJ) $(patsubst %.F90,$(ODIR)/%.o,$(SRC))
MODOBJ  = $(patsubst %.f90,$(ODIR)/%.o,$(SRCMOD)) 
MODOBJ := $(patsubst %.F90,$(ODIR)/%.o,$(MODOBJ))
MODMOD  = $(patsubst $(ODIR)/%.o,$(MDIR)/%.mod,$(MODOBJ))

dir_guard=@mkdir -p $(@D)

default: $(MODOBJ) $(OBJ)
	@mkdir -p $(BDIR)
	$(FC) $(MODOBJ) $(OBJ) -o $(BDIR)/$(PROGRAM) $(libpath)

#version.f90:
#	@./get_git_version.sh
	
# Compiler steps for all objects
$(ODIR)/%.o: %.F90
	@mkdir -p $(MDIR)
	$(dir_guard)
	$(FC) $(FCFLAGS) -o $@ $^ 

$(ODIR)/version.o: version.f90
	@mkdir -p $(MDIR)
	$(dir_guard)
	$(FC) $(FCFLAGS) -o $@ $^ 

$(ODIR)/%.o: %.f90
	@mkdir -p $(MDIR)
	$(dir_guard)
	$(FC) $(FCFLAGS) -o $@ $^ 

# Linker
$(BDIR)/$(PROGRAM) : $(OBJ)
	$(dir_guard)
	$(FC) $(FLFLAGS) -o $@ $< $(libpath) 

# Dependencies : Main program depends on all modules
$(ODIR)/$(PRGOBJ) : $(MODOBJ)

debug: 

	@echo "SRCF90  = $(SRCF90)"
	@echo "SRCMOD  = $(SRCMOD)"
	@echo "OBJ     = $(OBJ)"
	@echo "MODOBJ  = $(MODOBJ)"
	@echo "MODMOD  = $(MODMOD)"
	@echo "PROGRAM = $(PROGRAM)"
	@echo "PRGOBJ  = $(PRGOBJ)"
#	DEBUG=1 make clean $(BDIR)/$(PROGRAM)

install:
	cp $(BDIR)/$(PROGRAM) $(DIR)

clean:
	rm -f $(ODIR)/*.o $(MDIR)/*mod $(BDIR)/$(PROGRAM) $(BDIR)/debug

