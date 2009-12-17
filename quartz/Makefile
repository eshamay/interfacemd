include ../Makefile

%.o: %.cpp %.h
	$(CXX) $(CPPFLAGS) -c -o $@ $<

PDBSYSTEM		= $(MDSYSTEM) ../pdbfile.o

QUARTZ			= $(PDBSYSTEM) make-quartz-slab.cpp
MAKEONEMOL		= $(PDBSYSTEM) make-one-mol.cpp
CONNECTIVITY	= $(PDBSYSTEM) quartz-connectivity.cpp
TEMPFILES		= $(WATERSYSTEM) ../pdbfile.o temp.cpp

quartz: $(PDBSYSTEM) make-quartz-slab.cpp
	$(CXX) $(CPPFLAGS) $(PDBSYSTEM) make-quartz-slab.cpp -llapack -o make-quartz-slab

temp: $(TEMPFILES)
	$(CXX) $(CPPFLAGS) $(TEMPFILES) -llapack -o temp 

makeonemolecule: $(MAKEONEMOL)
	$(CXX) $(CPPFLAGS) $(MAKEONEMOL) -o make-one-mol 

connectivity: $(CONNECTIVITY)
	$(CXX) $(CPPFLAGS) $(CONNECTIVITY) -o connect

findoh: $(PDBSYSTEM) find-oh.cpp
	$(CXX) $(CPPFLAGS) $(PDBSYSTEM) find-oh.cpp -o find-oh

independent-mols: $(PDBSYSTEM) make-independent-mols.cpp
	$(CXX) $(CPPFLAGS) $(PDBSYSTEM) make-independent-mols.cpp -o make-ind-mols
	
findtypes: $(PDBSYSTEM) find-atom-types.cpp
	$(CXX) $(CPPFLAGS) $(PDBSYSTEM) find-atom-types.cpp -o find-types


  	

cleanquartz :
	rm -f *.o make-quartz-slab temp make-one-mol connect find-oh
