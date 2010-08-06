include ../Makefile

SO2ANGLEBONDLENGTH		= $(ANALYSIS) $(WATER) $(IONS) $(ORGANIC) so2-angle-bond-histograms.o
H2OANGLEBONDLENGTH		= $(ANALYSIS) $(WATER) $(IONS) $(ORGANIC) h2o-angle-bond-histograms.o

test : $(TEST) 
	icpc $(TEST) $(LIBS) $(LAPACKLIBS) -o test

so2-anglebondhisto : $(SO2ANGLEBONDLENGTH) 
	icpc $(SO2ANGLEBONDLENGTH) $(LIBS) -o ../bin/so2-angle-bond-histograms

h2o-anglebondhisto : $(H2OANGLEBONDLENGTH) 
	icpc $(H2OANGLEBONDLENGTH) $(LIBS) -o ../bin/h2o-angle-bond-histograms

cleanstructure :
	rm -f *.o 
