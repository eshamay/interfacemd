include ../Makefile

STRUCTURE				= $(ANALYSIS) $(WATER) $(IONS) $(ORGANIC) 
SYSTEMANALYSES			= $(STRUCTURE) so2-structure-analysis.o

test : $(TEST) 
	icpc $(TEST) $(LIBS) $(LAPACKLIBS) -o test

so2-structure-analyses : $(SYSTEMANALYSES)
	icpc $(SYSTEMANALYSES) $(LIBS) -o ../bin/so2-structure-analyses


cleanstructure :
	rm -f *.o 
