include ../Makefile

STRUCTURE				= $(ANALYSIS) $(WATER) $(IONS) $(ORGANIC) 
ATOMICDENSITY			= $(STRUCTURE) atomic-density-analysis.o

test : $(TEST) 
	icpc $(TEST) $(LIBS) $(LAPACKLIBS) -o test

atomic-density-analysis : $(ATOMICDENSITY)
	icpc $(ATOMICDENSITY) $(LIBS) -o ../bin/atomic-density-analysis


cleanstructure :
	rm -f *.o 
