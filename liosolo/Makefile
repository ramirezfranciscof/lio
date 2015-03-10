######################################################################
ifeq ($(ifort),1)
  FC   = ifort
  LIBS+= -module ../lioamber/obj
else
  FC   = gfortran
  LIBS+= -I../lioamber/obj -J../lioamber/obj
endif

LIBS += -L../g2g -lg2g
LIBS += -L../lioamber -llio-g2g
FLAGS =

ifeq ($(magma),1)
        DEFINE   += -Dmagma
        LIBS     += -L$(MAGMAROOT)/lib -lmagma
endif

<<<<<<< HEAD
ifeq ($(debug),1)
        DBUG = -g -O0 -debug all $(DEFINE)
endif

liosolo: main.f
	ifort $(DBUG) -L ../lioamber/ -llio-g2g -L ../g2g/ -lg2g  -I../lioamber/ $(LIBS)  main.f -o liosolo
=======
######################################################################
liosolo: main.f
	$(FC) $< -o $@ $(FLAGS) $(LIBS)

>>>>>>> 7d06646cb1333889dc332b2c15ba48ac06d0dcde
clean:
	rm -rf liosolo
######################################################################
