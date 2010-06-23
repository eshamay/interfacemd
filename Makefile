SRCLIB	 	= $(HOME)/md/src

BOOST		= $(HOME)/src/boost_1_43_0
LAPACK		= $(HOME)/src/lapack-3.2.1
LAPACKLIBS	= -lmkl_lapack -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

XDRLIB		= $(XDRDIR)/lib
XDRINC		= $(XDRDIR)/include

LIBS		= -lconfig++
CINCLUDE	= -I$(SRCLIB) -I$(HOME)/share/include -I$(BOOST) -L$(HOME)/share/lib #-I$(XDRINC)
CPPFLAGS	= #-L$(XDRDIR)/lib
#CXX			= mpiCC -g
CXXDEBUG	= -O0 -g3 -ggdb -D_GLIBCXX_DEBUG -Wno-deprecated #-wd981,1599,1572,383
CXXOPTIMIZE = -O2 -finline-functions -finline -funroll-loops

CXXFLAGS    = $(CINCLUDE) -Wall -Drestrict= -ftemplate-depth-100 $(CXXDEBUG) 
CXX			= g++ $(CXXFLAGS)

MPICXX		= mpiCC -g -I$(MPI)/include

MATH	= $(SRCLIB)/vecr.o $(SRCLIB)/matrixr.o 

MDSYSTEM = $(MATH) $(SRCLIB)/atom.o $(SRCLIB)/molecule.o $(SRCLIB)/mdsystem.o $(SRCLIB)/bond.o $(SRCLIB)/graph.o

WATER	= $(SRCLIB)/h2o.o 

IONS	= $(SRCLIB)/hno3.o $(SRCLIB)/h3o.o $(SRCLIB)/oh.o 

ORGANIC = $(SRCLIB)/decane.o $(SRCLIB)/carbonchain.o 

XYZSYSTEM = $(MDSYSTEM) $(SRCLIB)/xyzsystem.o $(SRCLIB)/xyzfile.o $(SRCLIB)/wannier.o 

AMBERSYSTEM	= $(MDSYSTEM) $(SRCLIB)/ambersystem.o $(SRCLIB)/crdfile.o $(SRCLIB)/forcefile.o $(SRCLIB)/topfile.o

GMXSYS	= $(MDSYSTEM) $(SRCLIB)/trrfile.o $(SRCLIB)/grofile.o $(SRCLIB)/gmxsystem.o

clean:
	rm -f *.o bin/*

cleanall:
	( make clean )
	( cd test ; make cleantest )
	( cd mpi ; make clean )
	( cd quartz; make cleanquartz )
	( cd sfg; make cleansfg )
	( cd system-analysis; make cleansysan )

%.o: %.cpp %.h
	$(CXX) $(CPPFLAGS) -c -o $@ $<

