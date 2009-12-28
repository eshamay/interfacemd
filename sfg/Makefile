include ../Makefile

%.o: %.cpp %.h
	$(CXX) $(CPPFLAGS) -c -o $@ $<

MORITAFILES		= $(WATERSYSTEM) watersfg.o morita.o 
MPIMORITAFILES	= $(MORITAFILES) ../mpi/mpisys.o
MORITA2002		= $(WATERSYSTEM) moritasfg2002.o morita2002.o

morita : $(MORITAFILES)
	$(CXX) $(CPPFLAGS) $(MORITAFILES) -llapack -o ../bin/morita

mpimorita : $(MPIMORITAFILES)
	$(MPICXX) $(CPPFLAGS) $(MPILIBS) $(MPIMORITAFILES) -o ../bin/mpi-morita


cleansfg :
	rm -f *.o
