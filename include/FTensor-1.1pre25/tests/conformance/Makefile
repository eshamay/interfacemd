all: test_all
VPATH=T2:T2s:T3:T3as:T3dg:T4ddg

CXX = g++-3.2
CXXOPTIMIZE =
CXXFLAGS=-ftemplate-depth-100 -Drestrict= -DFTENSOR_DEBUG $(CXXOPTIMIZE)
# CXXFLAGS=-O3 -ftemplate-depth-100 -Drestrict=__restrict__

# CXX = /home/boo/foreign_compilers/usr/local/KAI/KCC.pu-4.0d-1/KCC_BASE/bin/KCC
# # CXXFLAGS=+K3 --restrict --no_exceptions --inline_auto_space_time=750000000000000 --inline_implicit_space_time=200000000000000 --inline_generated_space_time=40000000000000.0 --inline_auto_space_time=100000000000000.0 --max_pending_instantiations 100
# CXXFLAGS= -g --restrict --no_exceptions --max_pending_instantiations 100 -DFTENSOR_DEBUG
# # CXXFLAGS= --restrict --no_exceptions --max_pending_instantiations 100
# # CXXFLAGS= -Drestrict=

# CXX = icc
# CXXFLAGS= -restrict -O0 -g -DFTENSOR_DEBUG -lcprts -lcxa -lguide -limf -lirc -lircmt -lsvml

# CXX = pgCC
# CXXFLAGS = -Drestrict= -DFTENSOR_DEBUG

# CXX = /home/boo/Compilers/Intel/opt/intel/compiler60/ia32/bin/icc
# CXXFLAGS = -restrict -DFTENSOR_DEBUG


TENSOR0=test_T0.o
TENSOR1=test_T1.o

TENSOR2=test_T2.o test_T2I.o test_T2II.o test_T2III.o test_T2IV.o test_T2V.o test_T2VI.o test_T2VII.o test_T2VIII.o test_T2IX.o test_T2X.o test_T2XI.o test_T2XII.o test_T2XIII.o test_T2XIV.o test_T2XV.o test_T2XVI.o test_T2XVII.o test_T2XVIII.o test_T2XIX.o test_T2XX.o test_T2XXI.o test_T2XXII.o test_T2XXIII.o test_T2XXIV.o test_T2XXV.o test_T2XXVI.o test_T2XXVII.o test_T2XXVIII.o test_T2XXIX.o test_T2XXX.o test_T2XXXI.o test_T2XXXII.o test_T2XXXIII.o test_T2XXXIV.o test_T2XXXV.o test_T2XXXVI.o test_T2XXXVII.o test_T2XXXVIII.o
TENSOR2S=test_T2s.o test_T2sI.o test_T2sII.o test_T2sIII.o test_T2sIV.o test_T2sV.o test_T2sVI.o test_T2sVII.o test_T2sVIII.o test_T2sIX.o test_T2sX.o test_T2sXI.o test_T2sXII.o test_T2sXIII.o test_T2sXIV.o test_T2sXV.o test_T2sXVI.o test_T2sXVII.o test_T2sXVIII.o test_T2sXIX.o test_T2sXX.o test_T2sXXI.o test_T2sXXII.o test_T2sXXIII.o test_T2sXXIV.o test_T2sXXV.o test_T2sXXVI.o test_T2sXXVII.o test_T2sXXVIII.o test_T2sXXIX.o test_T2sXXX.o test_T2sXXXI.o test_T2sXXXII.o test_T2sXXXIII.o test_T2sXXXIV.o test_T2sXXXV.o test_T2sXXXVI.o test_T2sXXXVII.o test_T2sXXXVIII.o test_T2sXXXIX.o test_T2sXL.o test_T2sXLI.o test_T2sXLII.o test_T2sXLIII.o test_T2sXLIV.o

TENSOR3=test_T3.o test_T3II.o
TENSOR3AS=test_T3as.o test_T3asI.o test_T3asII.o test_T3asIII.o test_T3asIV.o test_T3asV.o test_T3asVI.o test_T3asVII.o test_T3asVIII.o test_T3asIX.o test_T3asX.o test_T3asXI.o test_T3asXII.o test_T3asXIII.o test_T3asXIV.o
TENSOR3CH=test_T3ch.o
TENSOR3DG=test_T3dg.o test_T3dgI.o test_T3dgII.o test_T3dgIII.o test_T3dgIV.o test_T3dgV.o test_T3dgVI.o test_T3dgVII.o test_T3dgVIII.o test_T3dgIX.o test_T3dgX.o test_T3dgXI.o test_T3dgXII.o test_T3dgXIII.o test_T3dgXIV.o test_T3dgXV.o test_T3dgXVI.o test_T3dgXVII.o test_T3dgXVIII.o test_T3dgXIX.o  test_T3dgXX.o test_T3dgXXI.o test_T3dgXXII.o test_T3dgXXIII.o test_T3dgXXIV.o test_T3dgXXV.o test_T3dgXXVI.o test_T3dgXXVII.o test_T3dgXXVIII.o test_T3dgXXIX.o test_T3dgXXX.o test_T3dgXXXI.o test_T3dgXXXII.o test_T3dgXXXIII.o test_T3dgXXXIV.o test_T3dgXXXV.o test_T3dgCII.o test_T3dgCIII.o
TENSOR4DDG=test_T4ddg.o test_T4ddgII.o test_T4ddgIII.o test_T4ddgIV.o test_T4ddgV.o

testing_files=$(TENSOR0) $(TENSOR1) $(TENSOR2) $(TENSOR2S) $(TENSOR3) $(TENSOR3AS) $(TENSOR3CH) $(TENSOR3DG) $(TENSOR4DDG) testcases.o

test_all: $(testing_files)
	$(CXX) $(CXXFLAGS) -o test_all $(testing_files)

compiler_files=$(TENSOR0) $(TENSOR1) $(TENSOR2) $(TENSOR2S) $(TENSOR3AS) $(TENSOR3CH) testcases_compiler.o

test_compiler: $(compiler_files)
	$(CXX) $(CXXFLAGS) -o test_compiler $(compiler_files)

clean:
	rm -f *.o test_all test_compiler
