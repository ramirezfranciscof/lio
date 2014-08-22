LIBS=-L../g2g -lg2g -lstdc++ \
     -L $(CUDA_HOME)/lib64/ -lcudart -lcublas \
     -L ../lioamber -llio-g2g

liosolo:
	ifort $(LIBS) -I../lioamber main.f -o  liosolo
clean:
	rm -rf liosolo
 

