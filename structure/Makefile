include ../Makefile

STRUCTURE				= $(ANALYSIS) $(WATER) $(IONS) $(ORGANIC) 
MORITA					= ../sfg/moritah2o.o ../sfg/amber-morita2002.o 
SYSTEMANALYSES			= $(STRUCTURE) $(MORITA) angle-bond-analysis.o neighbor-analysis.o structure-analyzer.o 



test : $(TEST) 
	icpc $(TEST) $(LIBS) $(LAPACKLIBS) -o test

structure-analyzer : $(SYSTEMANALYSES)
	icpc $(SYSTEMANALYSES) $(LIBS) -o ../bin/structure-analyzer


cleanstructure :
	rm -f *.o 
