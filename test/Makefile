include ../Makefile

MDSYSTEM = ../atom.o ../h2o.o ../hno3.o ../molecule.o ../vecr.o ../matrixr.o ../bond.o ../oh.o ../h3o.o
XYZSYSTEM = ../xyzsystem.o ../xyzfile.o ../wannier.o ../graph.o
AMBERSYSTEM	= ../ambersystem.o ../crdfile.o ../forcefile.o ../topfile.o

WATERSYSTEM	= $(MDSYSTEM) $(XYZSYSTEM) $(AMBERSYSTEM)
CARBONCHAIN	= ../carbonchain.o ../decane.o
MORITAFILES	= morita.o ../watersfg.o $(WATERSYSTEM)
MPIMORITAFILES	= $(MORITAFILES) ../mpi/mpisys.o
MORITA2002	= $(WATERSYSTEM) ../moritasfg2002.o morita2002.o

DENSITYFILES = $(WATERSYSTEM) densitytest.o
DIPOLETCF = $(XYZSYSTEM) dipoletcf.o ../hno3analysis.o
QPT = $(XYZSYSTEM) ../hno3analysis.o Qpt.o
MDCRD2PDB = $(AMBERSYSTEM) ../pdbfile.o mdcrd2pdb.o
FILLBOX = $(AMBERSYSTEM) ../pdbfile.o ../boxfiller.o
MPIFILES	=	$(AMBERSYSTEM) ../mpi/mpisys.o mpitest.o
ORDERPARAMS = $(WATERSYSTEM) orderparams.o
ORDERPARAMSHISTO = $(WATERSYSTEM) orderparams-histo.o
TEMPFILES	= $(MDSYSTEM) $(XYZSYSTEM) ../pdbfile.o temp.cpp
COORDINATION	= $(WATERSYSTEM) coordination.o
ORDERCOORDS		= $(WATERSYSTEM) ordercoords.o
HCORIENTATION	= $(WATERSYSTEM) $(CARBONCHAIN) carbonchainsystem.o

ROTATE			= $(MDSYSTEM) $(XYZSYSTEM) ../pdbfile.o molecule-rotate.cpp

hcorientation: $(HCORIENTATION)
	$(CXX) $(CPPFLAGS) $(HCORIENTATION) -llapack -o carbon-orientation

ambertest : $(WATERSYSTEM) ambertest.o
	$(CXX) $(CPPFLAGS) $(WATERSYSTEM) ambertest.o -llapack -o ambertest

densitytest : $(DENSITYFILES)
	$(CXX) $(CPPFLAGS) $(DENSITYFILES) -llapack -o densitytest

morita : $(MORITAFILES)
	$(CXX) $(CPPFLAGS) $(MORITAFILES) -llapack -o morita

mpimorita : $(MPIMORITAFILES)
	$(MPICXX) $(CPPFLAGS) $(MPILIBS) $(MPIMORITAFILES) -o mpi-morita

matrixtest : matrixtest.o ../vecr.o ../matrixr.o
	$(CXX) $(CPPFLAGS) matrixtest.o ../vecr.o ../matrixr.o -o matrixtest

toptest : toptest.o ../topfile.o
	$(CXX) $(CPPFLAGS) toptest.o ../topfile.o -o toptest

angledistro: angledistro.o $(AMBERSYSTEM)
	$(CXX) $(CPPFLAGS) angledistro.o $(AMBERSYSTEM) -o angledistro

orderparams:  $(ORDERPARAMS)
	$(CXX) $(CPPFLAGS) $(ORDERPARAMS) -llapack -o orderparams

orderparams-histo:  $(ORDERPARAMSHISTO)
	$(CXX) $(CPPFLAGS) $(ORDERPARAMSHISTO) -llapack -o orderparams-histo

dipoletcf : $(DIPOLETCF)
	$(CXX) $(CPPFLAGS) $(DIPOLETCF) -o dipoletcf

qpt : $(QPT)
	$(CXX) $(CPPFLAGS) -o qpt_test $(QPT)

closest : $(AMBERSYSTEM) closest.o ../pdbfile.o
	$(CXX) $(CPPFLAGS) closest.o ../pdbfile.o $(AMBERSYSTEM) -o closest

mdcrd2pdb : $(MDCRD2PDB)
	$(CXX) $(CPPFLAGS) $(MDCRD2PDB) -o mdcrd2pdb

fillbox : $(FILLBOX) fillbox.o
	$(CXX) $(CPPFLAGS) fillbox.o $(FILLBOX) -o fillbox

temp: $(TEMPFILES)
	$(CXX) $(CPPFLAGS) -llapack -g3 -O0 -o temp $(TEMPFILES)

rotate: $(ROTATE)
	$(CXX) $(CPPFLAGS) -llapack -g3 -O0 -o rotate $(ROTATE)

sandbox: $(MDSYSTEM) $(AMBERSYSTEM)
	$(CXX) $(CPPFLAGS) -g3 -O0 -llapack -o matrix-sandbox $(MDSYSTEM) $(AMBERSYSTEM) matrix-sandbox.cpp

mpitest: $(MPIFILES)
	mpiCC $(CPPFLAGS) $(MPILIBS) -o mpitest $(MPIFILES)

coordination: $(COORDINATION)
	$(CXX) $(CPPFLAGS) $(COORDINATION) -o coordination

ordercoords: $(ORDERCOORDS)
	$(CXX) $(CPPFLAGS) $(ORDERCOORDS) -o ordercoords

morita2002: $(MORITA2002)
	$(CXX) $(CPPFLAGS) $(MORITA2002) -g3 -llapack -o morita2002

cleantest :
	rm -f *.o morita matrixtest toptest ambertest densitytest dipoletcf angledistro qpt_test a.out mdcrd2pdb closest fillbox test temp mpitest orderparams ordercoords morita2002

#densitytest.o: densitytest.cpp densitytest.h
#toptest.o: toptest.cpp toptest.h
#matrixtest.o: matrixtest.cpp ../matrixr.cpp ../matrixr.h
#morita.o: morita.cpp morita.h
#ambertest.o: ambertest.cpp ambertest.h
#angletest.o: angledistro.cpp angledistro.h
#Qpt.o: Qpt.cpp Qpt.h
#closest.o: closest.cpp
#mdcrd2pdb.o: mdcrd2pdb.cpp
#fillbox.o: fillbox.cpp
