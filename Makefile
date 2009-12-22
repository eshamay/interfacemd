SRCLIB	 	= $(HOME)/work/src
MPILIBS		= $(MPI)/lib/libmpi.so $(MPI)/lib/libmpi_cxx.so
FTENSOR		= $(SRCLIB)/include/FTensor-1.1pre25
BOOST		= /common/src/boost_1_41_0
CINCLUDE	= -I$(SRCLIB) -I/usr/include -I$(BOOST)
CLIBS		= -L$(MKL) -lmkl_lapack -lmkl -lguide -lpthread
CPPFLAGS	= $(CINCLUDE)
#CXX			= mpiCC -g
CXXDEBUG	= -g3 -ggdb -Wall -D_GLIBCXX_DEBUG
CXXOPTIMIZE = -O2 -finline-functions -finline-limit-1000 -funroll-loops
#CXXFLAGS	= -ftemplate-depth-100 -Drestrict= $(CXXOPTIMIZE)
CXXFLAGS    = -ftemplate-depth-100 -Drestrict= $(CXXDEBUG) #-wd981,1599,1572,383
CXX			= g++ $(CXXFLAGS)
#CXX			= icpc $(CXXFLAGS)
MPICXX		= mpiCC -g -I$(MPI)/include

MDSYSTEM = $(SRCLIB)/vecr.o $(SRCLIB)/matrixr.o $(SRCLIB)/atom.o $(SRCLIB)/molecule.o $(SRCLIB)/mdsystem.o $(SRCLIB)/h2o.o $(SRCLIB)/hno3.o $(SRCLIB)/bond.o $(SRCLIB)/oh.o $(SRCLIB)/h3o.o $(SRCLIB)/decane.o $(SRCLIB)/carbonchain.o 
XYZSYSTEM = $(SRCLIB)/xyzsystem.o $(SRCLIB)/xyzfile.o $(SRCLIB)/wannier.o $(SRCLIB)/graph.o
AMBERSYSTEM	= $(SRCLIB)/ambersystem.o $(SRCLIB)/crdfile.o $(SRCLIB)/forcefile.o $(SRCLIB)/topfile.o
WATERSYSTEM	= $(MDSYSTEM) $(XYZSYSTEM) $(AMBERSYSTEM) 

clean:
	rm -f *.o

cleanall:
	( make clean )
	( cd test ; make cleantest )
	( cd mpi ; make clean )
	( cd quartz; make cleanquartz )
	( cd sfg; make cleansfg )

%.o: %.cpp %.h
	$(CXX) $(CPPFLAGS) -c -o $@ $<

