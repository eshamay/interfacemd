include ../Makefile

MORITAFILES		= $(WATERSYSTEM) watersfg.o morita.o 
TEST		= ../mdmath.o ../vecr.o ../matrixr.o test.o
MPIMORITAFILES		= $(MORITAFILES) ../mpi/mpisys.o
MORITA2002		= $(ANALYSIS) $(WATER) $(IONS) $(ORGANIC) morita2002.o

XYZSFG		= $(ANALYSIS) $(WATER) $(IONS) $(ORGANIC) xyzsfg.o

morita : $(MORITAFILES)
	$(CXX) $(CPPFLAGS) $(MORITAFILES) -llapack -o ../bin/morita

#morita2002 : $(MORITA2002)
#gfortran $(MORITA2002) $(LIBS) -llapack -lf77blas -lcblas -latlas -o ../bin/morita2002

morita2002 : $(MORITA2002)
	icpc $(MORITA2002) $(LAPACKLIBS) $(LIBS) -o ../bin/morita2002

xyzsfg : $(XYZSFG)
	icpc $(XYZSFG) $(LIBS) -lfftw3 -o ../bin/xyzsfg
	
mpimorita : $(MPIMORITAFILES)
	$(MPICXX) $(CPPFLAGS) $(MPILIBS) $(MPIMORITAFILES) -o ../bin/mpi-morita

test : $(TEST)
	icpc $(LAPACKLIBS) $(TEST) -o test

cleansfg :
	rm -f *.o test
