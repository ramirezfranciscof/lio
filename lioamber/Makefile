######################################################################
# LIBLIO MAKEFILE
######################################################################
# This is the main file directing the compilation for the 
# library LIO. It was designed to simplify the customization
# of every aspect of this process and thus it should not be
# modified under normal code tinkering. Instead, one should
# attempt first to modify the auxiliary .mk files (beware
# not to change these files' names though).
makefiles := Makefile

# This should contain the list of all objects and modules
# to be compiled and their internal dependencies/relations.
# It should also specify the locations of the source files
# (.f files) for these targets and include all .mk files
# containing information on how to compile modules.
makefiles += make_depend.mk

# This should contain the default flags for compilation
# as well as the specific flags for those targets that
# require them, and optional flags to be defined during
# compilation by the user.
makefiles += make_flags.mk

# This should contain the default flags for final linking
# as well as the optional linking flags to be defined
# during compilation by the user. This also includes all
# external libraries that are to be linked to LIO.
makefiles += make_links.mk

#
<<<<<<< HEAD
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# OLD PARTS
#--------------------------------------------------------------------#
# sources and objects
MOD_SRCS=$(wildcard mods/*.f)
# COMMON_SRCS=$(wildcard *.f)
GARCHA_SRCS=$(wildcard *.f)

MOD_OBJ=$(MOD_SRCS:%.f=%.o)
GARCHA_OBJ = #auxsubs.o
GARCHA_OBJ += $(MOD_OBJ) $(COMMON_SRCS:%.f=%.o) $(GARCHA_SRCS:%.f=%.o)

#--------------------------------------------------------------------#
# THE OLD WAY TO LIBS
#LIBS=-L/usr/lib -L/usr/lib64 -L../g2g -lg2g -lstdc++ \
      -L$(MKLROOT)/lib/intel64 -I$(MKLROOT)/include $(PROFILE) \
      -lmkl_lapack95_lp64 -Wl,--start-group -lmkl_intel_lp64 \
      -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread \
      -lm -Wl,-rpath='$$ORIGIN/' -Wl,-rpath='$$ORIGIN/../g2g'
#--------------------------------------------------------------------#


# flags
ifeq ($(print), 1)
        DEFINE += -DPRINT_MATRICES
endif

ifeq ($(profile),1)
  PROFILE  = -pg
else
  PROFILE  =
endif

ifeq ($(cpu_recompute),0)
  CXXFLAGS += -DCPU_RECOMPUTE=0
else
  CXXFLAGS += -DCPU_RECOMPUTE=1
endif

ifeq ($(full_double),1)
  CXXFLAGS += -DFULL_DOUBLE=1
  NVCCFLAGS += -DFULL_DOUBLE=1
else
  CXXFLAGS += -DFULL_DOUBLE=0
  NVCCFLAGS += -DFULL_DOUBLE=0
endif

DEFINE   += -Dpack -fPIC
CFLAGS   = -g -mp1 -ip -O3 -c -cpp  $(DEFINE) $(PROFILE)
CFLAGS2  = -g -O1 -c -cpp  $(DEFINE) $(PROFILE)
CFLAGS3  = -g -O3 -c -cpp -parallel   $(DEFINE) $(PROFILE)
CPPFLAGS = -P -traditional-cpp

ifeq ($(magma),1)
	DEFINE   += -Dmagma
	LIBS     += -L$(MAGMAROOT)/lib -lmagma
endif

#eliminada   -lmkl_solver_lp64

ifeq ($(debug),1)
        DBUG = -g -O0 -debug all
        CFLAGS = $(DBUG) -ip -c -cpp $(DEFINE)
endif

# directories
G2G_CFLAGS=-DG2G

ifeq ($(ultima_garcha),1)
	G2G_CFLAGS += -DULTIMA_GARCHA
endif

ifeq ($(convergencia_garcha),1)
	G2G_CFLAGS += -DINT3N_GARCHA
endif

ifeq ($(ultima_g2g),1)
	CFLAGS += -DULTIMA_G2G
	CFLAGS2 += -DULTIMA_G2G
endif

## Targets: make ##

all: dirs liblio-g2g.so
#	@echo '---------------------------'
#	@echo $(MKLROOT)
#	@echo '---------------------------'
=======
######################################################################

# This indicates the main target.
all : liblio-g2g.so

# This indicates where to store the objects and modules produced
# by the compilation process. Vpath is not reliable for finding
# these files; obj_path shoud be added to objects explicitly.
obj_path := obj
vpath %.o   $(obj_path)
vpath %.mod $(obj_path)

# This indicates which objects are to be compiled, in which
# order, and where to look for its source files (or other
# dependencies, like .mk files)
objects   =
src_paths =
include make_depend.mk
vpath %.f  $(src_paths)
vpath %.mk $(src_paths)

# This indicates which flags to use during the compilation
# of the different objects and which flags to use during
# the final linking of the library. It also contains the
# options to make with different compilers.
# LIBS is commented because for the time being LIO may rely
# on information set via environment variables.
FC =
FFLAGS =
LFLAGS =
#LIBS   =
include make_flags.mk
include make_links.mk
>>>>>>> 7d06646cb1333889dc332b2c15ba48ac06d0dcde

#
######################################################################
# GENERAL RULES : The following rules have been written to be
# general and make use of variables defined along this and the
# other make files. They should not be modified under normal
# circumstances; modifications of added files or change in the
# compilation flags should be done through these variables.
liblio-g2g.so: $(objects:%.o=$(obj_path)/%.o)
	$(FC) $(LFLAGS) $^ $(LIBS) -o $@

$(obj_path)/%.o   :  %.f $(makefiles) |  $(obj_path)
	$(FC) $(FFLAGS) -c $< -o $@ 

$(obj_path)/%.mod :  %.f $(makefiles) |  $(obj_path)
	$(FC) $(FFLAGS) -c $< -o $(@:%.mod=%.o)

$(obj_path) :
	mkdir -p $@

.PHONY: clean
clean:
	rm -rf liblio-g2g.so $(obj_path)

######################################################################
