include ../Makefile

STRUCTURE				= $(ANALYSIS) $(WATER) $(IONS) $(ORGANIC) 
SYSTEMANALYSES			= $(STRUCTURE) angle-bond-analysis.o neighbor-analysis.o structure-analyzer.o 

test : $(TEST) 
	icpc $(TEST) $(LIBS) $(LAPACKLIBS) -o test

structure-analyzer : $(SYSTEMANALYSES)
	icpc $(SYSTEMANALYSES) $(LIBS) -o ../bin/structure-analyzer


cleanstructure :
	rm -f *.o 
