# CXX = /home/boo/foreign_compilers/usr/local/KAI/KCC.pu-4.0d-1/KCC_BASE/bin/KCC
# #CXXFLAGS=-Drestrict= --max_pending_instantiations 1000
# CXXFLAGS=--restrict --no_exceptions --max_pending_instantiations 100
# # CXXFLAGS=+K3 --restrict --no_exceptions --inline_auto_space_time=750000000000000 --inline_implicit_space_time=200000000000000 --inline_generated_space_time=40000000000000.0 --inline_auto_space_time=100000000000000.0 --max_pending_instantiations 100

CXX = g++
CXXOPTIMIZE=-O3 -finline-functions -finline-limit-1000 -funroll-loops
CXXFLAGS=-ftemplate-depth-100 -Drestrict= $(CXXOPTIMIZE)
# CXXFLAGS=-O3 -ftemplate-depth-100 -Drestrict= -fomit-frame-pointer -finline-functions -finline-limit-1000000000000 -ffast-math -fno-rtti -fno-exceptions
# CXXFLAGS=-ftemplate-depth-100 -Drestrict=__restrict__
# CXXFLAGS=-ansi -ftemplate-depth-100 -Drestrict=
# CXXFLAGS=-O3 -ftemplate-depth-100 -Drestrict=__restrict__


# CXX = icc
# # CXXFLAGS= -restrict -O0
# # CXXFLAGS= -restrict -O0 -Qoption,cpp,--pending_instantiations,1000
# # CXXFLAGS= -restrict -O3
# # CXXFLAGS= -restrict -O3 -tpp6 -xi
# # CXXFLAGS= -restrict -O3 -tpp6 -xi -ip
# CXXFLAGS= -restrict -O3 -tpp6 -xi -ipo
# #CXXFLAGS= -restrict -O0 -g -I/home/boo/foreign_compilers/opt/intel/compiler50/ia32/include -I/usr/local/lib/gcc-lib/i686-pc-linux-gnu/2.95.3/include
# #CXXFLAGS= -restrict -O0 -g -L/home/boo/foreign_compilers/opt/intel/compiler50/ia32/lib -lcprts
# #CXXFLAGS= -restrict -O0 -g -lcprts
# # CXXFLAGS= -restrict -I/home/boo/intel_compiler/opt/intel/compiler50/ia32/include -I/usr/local/lib/gcc-lib/i686-pc-linux-gnu/2.95.3/include -O3 -tpp6 -xi
# # # CXXFLAGS= -Drestrict= -I/home/boo/intel_compiler/opt/intel/compiler50/ia32/include -I/usr/lib/gcc-lib/i386-linux/2.95.2/include -O3 -tpp6 -xi -ip

# CXX = pgCC
# # CXXFLAGS = -Drestrict= -O2 -Munroll -Mnoframe  -Minline=levels:10 --no_exceptions -Mvect -Mvect=assoc -Mneginfo=loop

# CXXFLAGS = -Drestrict= --pending_instantiations=1000
# # CXXFLAGS = -Drestrict= -fast -Minline=levels:10 --no_exceptions

# # CXXFLAGS = -Drestrict= -O2 -Munroll -Mnoframe  -Minline=levels:100 --no_exceptions


% :: %,v

all: little littlefast example examplefast speed_test one_over one_over_fast


one_over: one_over_1_minus_x1 one_over_1_minus_x2 one_over_1_minus_x3 one_over_1_minus_x4 one_over_1_minus_x5 one_over_1_minus_x6 one_over_1_minus_x7 one_over_1_minus_x8 one_over_1_minus_x9
one_over_fast: one_over_1_minus_x_fast1 one_over_1_minus_x_fast2 one_over_1_minus_x_fast3 one_over_1_minus_x_fast4 one_over_1_minus_x_fast5 one_over_1_minus_x_fast6 one_over_1_minus_x_fast7 one_over_1_minus_x_fast8 one_over_1_minus_x_fast9

algebra3: algebra3_test1 algebra3_test2 algebra3_test3 algebra3_test4 algebra3_test5 algebra3_test6 algebra3_test7 algebra3_test8 algebra3_test9


clean:
	rm -f *.o little littlefast example examplefast speed_test one_over_1_minus_x1 one_over_1_minus_x2 one_over_1_minus_x3 one_over_1_minus_x4 one_over_1_minus_x5 one_over_1_minus_x6 one_over_1_minus_x7 one_over_1_minus_x8 one_over_1_minus_x9 one_over_1_minus_x_fast1 one_over_1_minus_x_fast2 one_over_1_minus_x_fast3 one_over_1_minus_x_fast4 one_over_1_minus_x_fast5 one_over_1_minus_x_fast6 one_over_1_minus_x_fast7 one_over_1_minus_x_fast8 one_over_1_minus_x_fast9 algebra3_test1 algebra3_test2 algebra3_test3 algebra3_test4 algebra3_test5 algebra3_test6 algebra3_test7 algebra3_test8 algebra3_test9
