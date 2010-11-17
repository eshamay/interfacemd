include ../Makefile

PDBSYSTEM		= $(MDSYSTEM) ../pdbsystem.o

QUARTZ			= $(PDBSYSTEM) ../xyzfile.o quartz-slab.o
ROTATEUNITCELL			= $(PDBSYSTEM) ../xyzfile.o rotate-unitcell.o
MAKEONEMOL		= $(PDBSYSTEM) make-one-mol.cpp
CONNECTIVITY		= $(PDBSYSTEM) connectivity.cpp
TEMPFILES		= $(WATERSYSTEM) ../pdbfile.o temp.cpp
MOLGRID			= $(PDBSYSTEM) molgrid.cpp

quartz: $(QUARTZ)
	$(CXX) $(QUARTZ) -o ../bin/make-quartz-slab

rotateunitcell: $(ROTATEUNITCELL)
	$(CXX) $(ROTATEUNITCELL) -o ../bin/rotateunitcell

molgrid: $(MOLGRID)
	$(CXX) $(CPPFLAGS) $(MOLGRID) -llapack -o ../bin/molgrid

temp: $(TEMPFILES)
	$(CXX) $(CPPFLAGS) $(TEMPFILES) -llapack -o ../bin/temp 

makeonemolecule: $(MAKEONEMOL)
	$(CXX) $(CPPFLAGS) $(MAKEONEMOL) -o ../bin/make-one-mol 

connectivity: $(CONNECTIVITY)
	$(CXX) $(CPPFLAGS) $(CONNECTIVITY) -llapack -o ../bin/connect

findoh: $(PDBSYSTEM) find-oh.cpp
	$(CXX) $(CPPFLAGS) $(PDBSYSTEM) find-oh.cpp -o ../bin/find-oh

independent-mols: $(PDBSYSTEM) make-independent-mols.cpp
	$(CXX) $(CPPFLAGS) $(PDBSYSTEM) make-independent-mols.cpp -o ../bin/make-ind-mols
	
findtypes: $(PDBSYSTEM) find-atom-types.cpp
	$(CXX) $(CPPFLAGS) $(PDBSYSTEM) find-atom-types.cpp -llapack -o ../bin/find-types

cleanquartz :
	rm -f *.o ../bin/find-types ../bin/make-ind-mols a.out ../bin/connect ../bin/molgrid
