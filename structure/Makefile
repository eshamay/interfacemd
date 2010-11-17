include ../Makefile

STRUCTURE				= $(ANALYSIS) $(WATER) $(IONS) $(ORGANIC) 
MORITA					= ../sfg/moritah2o.o ../sfg/amber-morita2002.o 
SYSTEMANALYSES			= $(STRUCTURE) neighbor-analysis.o structure-analyzer.o #$(MORITA) 


test : $(TEST) 
	icpc $(TEST) $(LIBS) $(LAPACKLIBS) -o test

structure-analyzer : $(SYSTEMANALYSES)
	$(CXX) $(SYSTEMANALYSES) $(LAPACKLIBS) $(LIBS) -lxdrfile -o ../bin/structure-analyzer

cleanstructure :
	rm -f *.o 
