SRCLIB	 	= $(HOME)/md/src

CXX			= icpc -wd981,444,383,177

DEBUG		= -O0 -g3 -ggdb -D_GLIBCXX_DEBUG -Wno-deprecated -debug #-wd981,1599,1572,383
OPTIMIZE 	= -fast -finline-functions -finline -funroll-loops -m64
#CPPFLAGS    = -Wall -Drestrict= -ftemplate-depth-100 $(DEBUG) -L$(HOME)/share/lib
CPPFLAGS    = -Wall -ftemplate-depth-100 $(DEBUG)

LIBS		= -L$(HOME)/share/lib -L$(MPI_HOME)/lib -L$(ATLAS)/lib -lconfig++

ATLAS		= $(HOME)/share/atlas
BOOST		= $(HOME)/src/boost_1_43_0
LAPACK		= $(HOME)/src/lapack-3.2.1
LAPACKLIBS	= -lmkl_lapack -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
SCALAPACK	= -openmp -Wl,--start-group -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_intel_thread -Wl,--end-group -lpthread -lmpi -lm

XDRLIB		= $(XDRDIR)/lib
XDRINC		= $(XDRDIR)/include

#CPATH		= :$(SRCLIB):$(HOME)/share/include:$(ATLAS)/include:$(BOOST)#:$(XDRINC)


MATH	= $(SRCLIB)/vecr.o $(SRCLIB)/matrixr.o 

MDSYSTEM = $(MATH) $(SRCLIB)/atom.o $(SRCLIB)/molecule.o $(SRCLIB)/mdsystem.o

WATER	= $(SRCLIB)/h2o.o 

IONS	= $(SRCLIB)/hno3.o $(SRCLIB)/h3o.o $(SRCLIB)/oh.o $(SRCLIB)/so2.o

ORGANIC = $(SRCLIB)/carbonchain.o $(SRCLIB)/decane.o 

XYZSYSTEM = $(SRCLIB)/xyzsystem.o $(SRCLIB)/xyzfile.o $(SRCLIB)/wannier.o 

AMBERSYSTEM	= $(SRCLIB)/ambersystem.o $(SRCLIB)/crdfile.o $(SRCLIB)/forcefile.o $(SRCLIB)/topfile.o

GMXSYS	= $(SRCLIB)/trrfile.o $(SRCLIB)/grofile.o $(SRCLIB)/gmxsystem.o

WATERSYSTEM = $(MDSYSTEM) $(AMBERSYSTEM) $(XYZSYSTEM) $(SRCLIB)/graph.o 

ANALYSIS = $(WATERSYSTEM)

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

