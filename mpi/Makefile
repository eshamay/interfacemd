include ../Makefile

AMBERFILES	= ../vecr.o ../atom.o ../ambersystem.o ../connectmatrix.o ../crdfile.o ../dipoleparm.o ../forcefile.o ../h2o.o ../hno3.o ../molecule.o ../topfile.o
MPISYS	=	$(AMBERSYSTEM) mpisys.o

clean:
	rm -f *.o

%.o: %.cpp
	$(CPP) $(CPPFLAGS) -c -o $@ $<
