include ../Makefile

MORITAFILES		= $(WATERSYSTEM) watersfg.o morita.o 
TEST		= ../vecr.o ../matrixr.o test.o
MPIMORITAFILES		= $(MORITAFILES) ../mpi/mpisys.o
MORITA2002		= $(AMBERSYSTEM) $(WATER) $(IONS) $(ORGANIC) morita2002.o

morita : $(MORITAFILES)
	$(CXX) $(CPPFLAGS) $(MORITAFILES) -llapack -o ../bin/morita

morita2002 : $(MORITA2002)
	$(CXX) $(MORITA2002) $(LIBS) -o ../bin/morita2002

mpimorita : $(MPIMORITAFILES)
	$(MPICXX) $(CPPFLAGS) $(MPILIBS) $(MPIMORITAFILES) -o ../bin/mpi-morita

test : $(TEST)
	$(CXX) $(CPPFLAGS) $(TEST) -o test

cleansfg :
	rm -f *.o test
