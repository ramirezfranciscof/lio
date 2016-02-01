######################################################################
# LIOSOLO MAKEFILE
######################################################################
LIBS=
FLAGS = -fopenmp

ifeq ($(intel),1)
FLAGS += -fpp
else
FLAGS += -cpp
endif

ifeq ($(xeon_phi),1)
  FLAGS += -mmic
  LIBS= -Qoption,link,-rpath=/opt/intel/composer_xe_2013_sp1.3.174/mkl/lib/mic -Qoption,link,-rpath=/opt/intel/composer_xe_2013_sp1.3.174/compiler/lib/mic/
endif

ifeq ($(intel),1)
  FC   = ifort
  LIBS+= -module ../lioamber/obj -I../lioamber/obj 
else
  FC   = gfortran
  LIBS+= -I../lioamber/obj -J../lioamber/obj
endif

LIBS += -L../g2g -lg2g
LIBS += -L../lioamber -llio-g2g

ifeq ($(magma),1)
        DEFINE   += -Dmagma
        LIBS     += -L$(MAGMAROOT)/lib -lmagma
endif


######################################################################
all: liosolo liomd.x

liosolo: main.f
	$(FC) $< -o $@ $(FLAGS) $(LIBS)

liomd.x: liomd.f90
	$(FC) $< -o $@ $(FLAGS) $(LIBS)

clean:
	rm -rf liosolo

######################################################################
