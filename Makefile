#MDSRC	 	= $(HOME)/md/src


CXX			= icpc  -I$(EIGEN) -I$(BOOST) -I.. -wd981,444,383,177,1418,1782,869

DEBUG		= -O0 -g3 -ggdb -D_GLIBCXX_DEBUG -Wno-deprecated -DNDEBUG -debug #-wd981,1599,1572,383
OPTIMIZE 	= -finline-functions -finline -funroll-all-loops -O3 -DNDEBUG -m64 -fast -restrict
#CPPFLAGS    = -Wall -Drestrict= -ftemplate-depth-100 $(DEBUG) -L$(HOME)/share/lib
CPPFLAGS    = -Wall -ftemplate-depth-100 $(OPTIMIZE)

LIBS		= -L$(HOME)/share/lib -lconfig++ #-L$(MPI_HOME)/lib -L$(ATLAS)/lib 

LAPACKLIBS	= -lmkl_lapack -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
LAPACKLIBS32 = -lmkl_lapack -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
SCALAPACK	= -openmp -Wl,--start-group -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_intel_thread -Wl,--end-group -lpthread -lmpi -lm

XDRLIB		= $(XDRDIR)/lib
XDRINC		= $(XDRDIR)/include

#CPATH		= :$(MDSRC):$(HOME)/share/include:$(ATLAS)/include:$(BOOST)#:$(XDRINC)


MDSYSTEM = $(MDSRC)/atom.o $(MDSRC)/molecule.o $(MDSRC)/mdsystem.o

WATER	= $(MDSRC)/h2o.o 

IONS	= $(MDSRC)/hno3.o $(MDSRC)/h3o.o $(MDSRC)/oh.o $(MDSRC)/so2.o $(MDSRC)/h.o

ORGANIC = $(MDSRC)/alkane.o $(MDSRC)/decane.o 

XYZSYSTEM = $(MDSRC)/xyzsystem.o $(MDSRC)/xyzfile.o $(MDSRC)/wannier.o 

AMBERSYSTEM	= $(MDSRC)/ambersystem.o $(MDSRC)/crdfile.o $(MDSRC)/forcefile.o $(MDSRC)/topfile.o

GMXSYS	= $(MDSRC)/trrfile.o $(MDSRC)/grofile.o $(MDSRC)/gmxsystem.o

WATERSYSTEM = $(MDSYSTEM) $(AMBERSYSTEM) $(XYZSYSTEM) $(MDSRC)/graph.o

ANALYSIS = $(WATERSYSTEM) $(MDSRC)/dataoutput.o

clean:
	rm -f *.o

cleanall:
	( make clean )
	( cd test ; make cleantest )
	( cd mpi ; make clean )
	( cd quartz; make cleanquartz )
	( cd sfg; make cleansfg )
	( cd system-analysis; make cleansysan )

%.o: %.cpp %.h
	$(CXX) $(CPPFLAGS) -c -o $@ $<

