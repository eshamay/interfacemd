include ../Makefile

MORITAFILES		= $(WATERSYSTEM) watersfg.o morita.o 
TEST		= ../mdmath.o ../vecr.o ../matrixr.o test.o
MPIMORITAFILES		= $(MORITAFILES) ../mpi/mpisys.o
MORITA2002		= $(ANALYSIS) $(WATER) $(IONS) $(ORGANIC) moritah2o.o
FOO			= $(MDSYSTEM) $(WATER) foo.o
XYZIR		= $(ANALYSIS) $(WATER) $(IONS) $(ORGANIC) xyz_ir.o

morita : $(MORITAFILES)
	$(CXX) $(CPPFLAGS) $(MORITAFILES) -llapack -o ../bin/morita

#morita2002 : $(MORITA2002)
#gfortran $(MORITA2002) $(LIBS) -llapack -lf77blas -lcblas -latlas -o ../bin/morita2002

foo : $(FOO)
	icpc $(FOO) $(LAPACKLIBS) $(LIBS) -o foo

amber-morita2002 : $(MORITA2002) amber-morita2002.o
	icpc $(MORITA2002) amber-morita2002.o $(LAPACKLIBS) $(LIBS) -o ../bin/amber-morita2002 

cp2k-morita2002 : $(MORITA2002) cp2k-morita2002.o
	icpc $(MORITA2002) cp2k-morita2002.o $(LAPACKLIBS) $(LIBS) -o ../bin/cp2k-morita2002 

xyzir : $(XYZIR)
	icpc $(XYZIR) $(LIBS) -o ../bin/xyz_so2_ir
	
mpimorita : $(MPIMORITAFILES)
	$(MPICXX) $(CPPFLAGS) $(MPILIBS) $(MPIMORITAFILES) -o ../bin/mpi-morita

test : $(TEST)
	icpc $(LAPACKLIBS) $(TEST) -o test

cleansfg :
	rm -f *.o test
