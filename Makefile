SRCLIB	 	= $(HOME)/md/src

BOOST		= $(HOME)/src/boost-1_43_0
LAPACK		= $(HOME)/src/lapack-3.2.1
LAPACKLIBS	= -lmkl_lapack -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

XDRLIB		= $(XDRDIR)/lib
XDRINC		= $(XDRDIR)/include

CINCLUDE	= -I$(SRCLIB) -I$(HOME)/share/include -I/usr/include -I/usr/local/include -I$(BOOST) #-I$(XDRINC)
CPPFLAGS	= $(CINCLUDE) -L$(HOME)/share/lib -L$(LAPACK) -lconfig++ #-L$(XDRDIR)/lib
#CXX			= mpiCC -g
CXXDEBUG	= -g3 -ggdb -Wall -D_GLIBCXX_DEBUG
CXXOPTIMIZE = -O2 -finline-functions -finline -funroll-loops
#CXXFLAGS	= -ftemplate-depth-100 -Wall -Drestrict= $(CXXOPTIMIZE)
CXXFLAGS    = -ftemplate-depth-100 -Drestrict= -Wall -Wno-deprecated $(CXXDEBUG) #-wd981,1599,1572,383
CXX			= g++ $(CXXFLAGS)
MPICXX		= mpiCC -g -I$(MPI)/include

MDSYSTEM = $(SRCLIB)/vecr.o $(SRCLIB)/tensor.o $(SRCLIB)/matrixr.o $(SRCLIB)/atom.o $(SRCLIB)/molecule.o $(SRCLIB)/mdsystem.o $(SRCLIB)/h2o.o $(SRCLIB)/hno3.o $(SRCLIB)/bond.o $(SRCLIB)/oh.o $(SRCLIB)/h3o.o $(SRCLIB)/decane.o $(SRCLIB)/carbonchain.o 
XYZSYSTEM = $(SRCLIB)/xyzsystem.o $(SRCLIB)/xyzfile.o $(SRCLIB)/wannier.o $(SRCLIB)/graph.o
AMBERSYSTEM	= $(SRCLIB)/ambersystem.o $(SRCLIB)/crdfile.o $(SRCLIB)/forcefile.o $(SRCLIB)/topfile.o
GMXSYS	= $(SRCLIB)/trrfile.o $(SRCLIB)/grofile.o $(SRCLIB)/gmxsystem.o
WATERSYSTEM	= $(MDSYSTEM) $(XYZSYSTEM) $(AMBERSYSTEM)

clean:
	rm -f *.o bin/*

cleanall:
	( make clean )
	( cd test ; make cleantest )
	( cd mpi ; make clean )
	( cd quartz; make cleanquartz )
	( cd sfg; make cleansfg )
	( cd system-analysis; make cleansysan )
	( rm bin/* )

%.o: %.cpp %.h
	$(CXX) $(CPPFLAGS) -c -o $@ $<

