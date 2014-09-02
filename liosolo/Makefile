ifeq ($(magma),1)
        DEFINE   += -Dmagma
        LIBS     += -L$(MAGMAROOT)/lib -lmagma
endif

ifeq ($(debug),1)
        DBUG = -g -O0 -debug all $(DEFINE)
endif

liosolo: main.f
	ifort $(DBUG) -L ../lioamber/ -llio-g2g -L ../g2g/ -lg2g  -I../lioamber/ $(LIBS)  main.f -o liosolo
clean:
	rm liosolo

