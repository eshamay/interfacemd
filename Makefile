SRCLIB	 	= $(HOME)/work/src
MPILIBS		= $(MPI)/lib/libmpi.so $(MPI)/lib/libmpi_cxx.so 
CINCLUDE	= -I$(SRCLIB)
CLIBS		= -L$(MKL) -lmkl_lapack -lmkl -lguide -lpthread
CPPFLAGS	= $(CINCLUDE)
#CXX			= mpiCC -g
CXX			= g++ -g
MPICXX		= mpiCC -g -I$(MPI)/include

ANALYSISFILES	=	analysis.o xyzsystem.o connectmatrix.o dipoleparm.o h2o.o hno3.o matrixr.o molecule.o atom.o vecr.o wannier.o xyzfile.o

analysis: $(ANALYSISFILES)
	$(CXX) $(CFLAGS) $(ANALYSISFILES) -o analysis

clean:
	rm -f *.o a.out

cleanall:
	( make clean )
	( cd test ; make cleantest )
	( cd mpi ; make clean )

%.o: %.cpp %.h
	$(CXX) $(CPPFLAGS) -c -o $@ $<

#molecule.o: molecule.cpp molecule.h
#forcefile.o: forcefile.cpp forcefile.h
#crdfile.o: crdfile.cpp crdfile.h
#xyzsystem.o: xyzsystem.cpp xyzsystem.h
#pdbfile.o: pdbfile.cpp pdbfile.h atom.cpp atom.h molecule.cpp molecule.h
#topfile.o: topfile.cpp topfile.h
#density.o: density.cpp density.h
#matrixr.o: matrixr.cpp matrixr.h
#vecr.o: vecr.cpp vecr.h
#xyzfile.o: xyzfile.cpp xyzfile.h
#wannier.o: wannier.cpp wannier.h
#connectmatrix.o: connectmatrix.cpp connectmatrix.h
#watersfg.o: watersfg.cpp watersfg.h
#analysis.o: analysis.cpp
#dipolefieldtensor.o: dipolefieldtensor.cpp dipolefieldtensor.h
#atom.o: atom.cpp atom.h
#h2o.o: h2o.cpp h2o.h
#ambersystem.o: ambersystem.cpp ambersystem.h crdfile.o forcefile.o conntectmatrix.o topfile.o atom.o molecule.o hno3.o h2o.o
#dipoleparm.o: dipoleparm.cpp dipoleparm.h
#hno3.o: hno3.cpp hno3.h
#histogram.o: histogram.cpp histogram.h
#hno3analysis.o: hno3analysis.cpp hno3analysis.h
#boxfiller.o: boxfiller.cpp boxfiller.h
