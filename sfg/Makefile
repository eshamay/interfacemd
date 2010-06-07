include ../Makefile

MORITAFILES		= $(WATERSYSTEM) watersfg.o morita.o 
TEST		= $(WATERSYSTEM) test.o
MPIMORITAFILES		= $(MORITAFILES) ../mpi/mpisys.o
MORITA2002		= $(WATERSYSTEM) moritasfg2002.o morita2002.o

morita : $(MORITAFILES)
	$(CXX) $(CPPFLAGS) $(MORITAFILES) -llapack -lxdrfile -o ../bin/morita

mpimorita : $(MPIMORITAFILES)
	$(MPICXX) $(CPPFLAGS) $(MPILIBS) $(MPIMORITAFILES) -o ../bin/mpi-morita

test : $(TEST)
	$(CXX) $(CPPFLAGS) $(TEST) -llapack -lxdrfile -o test

cleansfg :
	rm -f *.o
