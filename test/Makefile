include ../Makefile

AMBERFILES	= ../vecr.o ../atom.o ../ambersystem.o ../crdfile.o ../forcefile.o ../h2o.o ../hno3.o ../molecule.o ../topfile.o ../matrixr.o
WATERSYSTEM	= $(AMBERFILES) ../watersystem.o ../adjacencymatrix.o ../bond.o

MORITAFILES	= morita.o ../watersfg.o $(WATERSYSTEM)
MPIMORITAFILES	= $(MORITAFILES) ../mpi/mpisys.o
DENSITYFILES = $(WATERSYSTEM) densitytest.o
XYZSYSTEM = ../xyzsystem.o ../xyzfile.o ../atom.o ../h2o.o ../hno3.o ../molecule.o ../vecr.o ../wannier.o ../bondgraph.o ../matrixr.o
DIPOLETCF = $(XYZSYSTEM) dipoletcf.o ../hno3analysis.o
QPT = $(XYZSYSTEM) ../hno3analysis.o Qpt.o
MDCRD2PDB = $(AMBERFILES) ../pdbfile.o mdcrd2pdb.o
FILLBOX = $(AMBERFILES) ../pdbfile.o ../boxfiller.o
MPIFILES	=	$(AMBERFILES) ../mpi/mpisys.o mpitest.o
ORDERPARAMS = $(WATERSYSTEM) orderparams.o
TEMPFILES	= temp.cpp $(XYZSYSTEM) ../moritasfg2002.o
COORDINATION	= $(WATERSYSTEM) coordination.o
ORDERCOORDS		= $(WATERSYSTEM) ordercoords.o

ambertest : ambertest.o $(AMBERFILES)
	$(CXX) $(CPPFLAGS) ambertest.o $(AMBERFILES) -o ambertest

densitytest : $(DENSITYFILES)
	$(CXX) $(CPPFLAGS) $(DENSITYFILES) -o densitytest

morita : $(MORITAFILES)
	$(CXX) $(CPPFLAGS) $(MORITAFILES) -o morita

mpimorita : $(MPIMORITAFILES)
	$(MPICXX) $(CPPFLAGS) $(MPILIBS) $(MPIMORITAFILES) -o mpi-morita

matrixtest : matrixtest.o ../vecr.o ../matrixr.o
	$(CXX) $(CPPFLAGS) matrixtest.o ../vecr.o ../matrixr.o -o matrixtest

toptest : toptest.o ../topfile.o
	$(CXX) $(CPPFLAGS) toptest.o ../topfile.o -o toptest

angledistro: angledistro.o $(AMBERFILES)
	$(CXX) $(CPPFLAGS) angledistro.o $(AMBERFILES) -o angledistro

orderparams:  $(ORDERPARAMS)
	$(CXX) $(CPPFLAGS) $(ORDERPARAMS) -o orderparams

dipoletcf : $(DIPOLETCF)
	$(CXX) $(CPPFLAGS) $(DIPOLETCF) -o dipoletcf

qpt : $(QPT)
	$(CXX) $(CPPFLAGS) -o qpt_test $(QPT)

closest : $(AMBERFILES) closest.o ../pdbfile.o
	$(CXX) $(CPPFLAGS) closest.o ../pdbfile.o $(AMBERFILES) -o closest

mdcrd2pdb : $(MDCRD2PDB)
	$(CXX) $(CPPFLAGS) $(MDCRD2PDB) -o mdcrd2pdb

fillbox : $(FILLBOX) fillbox.o
	$(CXX) $(CPPFLAGS) fillbox.o $(FILLBOX) -o fillbox

temp: $(TEMPFILES)
	$(CXX) $(TEMPFLAGS) -g -o temp $(TEMPFILES)

mpitest: $(MPIFILES)
	mpiCC $(CPPFLAGS) $(MPILIBS) -o mpitest $(MPIFILES)

coordination: $(COORDINATION)
	$(CXX) $(CPPFLAGS) $(COORDINATION) -o coordination

ordercoords: $(ORDERCOORDS)
	$(CXX) $(CPPFLAGS) $(ORDERCOORDS) -o ordercoords

cleantest :
	rm -f *.o morita matrixtest toptest ambertest densitytest dipoletcf angledistro qpt_test a.out mdcrd2pdb closest fillbox test temp mpitest orderparams ordercoords

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
