include ../Makefile

SYSTEM				= $(ANALYSIS) $(WATER) $(IONS) $(ORGANIC) 

MORITAFILES			= $(WATERSYSTEM) watersfg.o morita.o 
TEST				= ../mdmath.o ../vecr.o ../matrixr.o test.o
MPIMORITAFILES		= $(MORITAFILES) ../mpi/mpisys.o
MORITA2002			= $(SYSTEM) moritah2o.o
FOO					= $(MDSYSTEM) $(WATER) foo.o
XYZDIPOLETCF		= $(SYSTEM) xyz-system-dipole-tcf.o

LINKER = g++

morita : $(MORITAFILES)
	$(CXX) $(CPPFLAGS) $(MORITAFILES) -llapack -o ../bin/morita

foo : $(FOO)
	$(LINKER) $(FOO) $(LAPACKLIBS) $(LIBS) -o foo

amber-morita2002 : $(MORITA2002) amber-morita2002.o
	$(LINKER) $(MORITA2002) amber-morita2002.o $(LAPACKLIBS) $(LIBS) -o ../bin/amber-morita2002 

cp2k-morita2002 : $(MORITA2002) cp2k-morita2002.o
	$(LINKER) $(MORITA2002) cp2k-morita2002.o $(LAPACKLIBS32) $(LIBS) -o ../bin/cp2k-morita2002 

xyz-system-dipole-tcf : $(XYZDIPOLETCF)
	$(LINKER) $(XYZDIPOLETCF) $(LIBS) -o ../bin/xyz-system-dipole
	
mpimorita : $(MPIMORITAFILES)
	$(MPICXX) $(CPPFLAGS) $(MPILIBS) $(MPIMORITAFILES) -o ../bin/mpi-morita

test : $(TEST)
	$(LINKER) $(LAPACKLIBS) $(TEST) -o test

cleansfg :
	rm -f *.o test
