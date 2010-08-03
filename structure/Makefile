include ../Makefile

SO2ANGLEBONDLENGTH		= $(ANALYSIS) $(WATER) $(IONS) $(ORGANIC) so2-angle-bondlength.o
H2OANGLEBONDLENGTH		= $(ANALYSIS) $(WATER) $(IONS) $(ORGANIC) h2o-angle-bondlength.o

test : $(TEST) 
	icpc $(TEST) $(LIBS) $(LAPACKLIBS) -o test

so2-anglebond : $(SO2ANGLEBONDLENGTH) 
	icpc $(SO2ANGLEBONDLENGTH) $(LIBS) -o ../bin/so2-angle-bondlength

h2o-anglebond : $(H2OANGLEBONDLENGTH) 
	icpc $(H2OANGLEBONDLENGTH) $(LIBS) -o ../bin/h2o-angle-bondlength

cleanstructure :
	rm -f *.o 
