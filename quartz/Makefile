include ../Makefile

PDBSYSTEM		= $(MDSYSTEM) ../pdbfile.o

QUARTZ			= $(PDBSYSTEM) make-quartz-slab.cpp
MAKEONEMOL		= $(PDBSYSTEM) make-one-mol.cpp
CONNECTIVITY	= $(PDBSYSTEM) quartz-connectivity.cpp
TEMPFILES		= $(WATERSYSTEM) ../pdbfile.o temp.cpp

quartz: $(PDBSYSTEM) make-quartz-slab.cpp
	$(CXX) $(CPPFLAGS) $(PDBSYSTEM) make-quartz-slab.cpp -llapack -o ../bin/make-quartz-slab

temp: $(TEMPFILES)
	$(CXX) $(CPPFLAGS) $(TEMPFILES) -llapack -o ../bin/temp 

makeonemolecule: $(MAKEONEMOL)
	$(CXX) $(CPPFLAGS) $(MAKEONEMOL) -o ../bin/make-one-mol 

connectivity: $(CONNECTIVITY)
	$(CXX) $(CPPFLAGS) $(CONNECTIVITY) -o ../bin/connect

findoh: $(PDBSYSTEM) find-oh.cpp
	$(CXX) $(CPPFLAGS) $(PDBSYSTEM) find-oh.cpp -o ../bin/find-oh

independent-mols: $(PDBSYSTEM) make-independent-mols.cpp
	$(CXX) $(CPPFLAGS) $(PDBSYSTEM) make-independent-mols.cpp -o ../bin/make-ind-mols
	
findtypes: $(PDBSYSTEM) find-atom-types.cpp
	$(CXX) $(CPPFLAGS) $(PDBSYSTEM) find-atom-types.cpp -llapack -o ../bin/find-types


  	

cleanquartz :
	rm -f *.o
