#*********************************************************************************************
__MAC__ ?= false
#*********************************************************************************************
ifeq ($(__MAC__),false)
	CXX        = g++
	OMPFLAG    = -fopenmp
	CXXFLAGS   = -fopenmp -Wall -std=c++20 -O2#  -fsanitize=address
	MPFRFLAGS  = -lmpfr -lgmp
else
	CXX        = clang++
	CXXFLAGS   = -arch arm64 -std=c++20 -march=native -funroll-loops -I/opt/homebrew/Cellar/ -I/usr/local/include -I/opt/homebrew/Cellar/mpfr/4.2.2/include/ -O2 -D __MAC__
	MPFRFLAGS  = -L/usr/local/lib -lmpfr -lgmp
endif

#*********************************************************************************************

all: def_demo grid_demo local_demo global_demo

#*********************************************************************************************

def_demo: deflection_demo.o functions.o deflection.o 
ifeq ($(__MAC__),false)
	$(CXX) $(CXXFLAGS) -o def_demo deflection_demo.o functions.o deflection.o $(MPFRFLAGS)
else
	$(CXX) $(CXXFLAGS) -o def_demo deflection_demo.o functions.o deflection.o $(MPFRFLAGS)
endif

#*********************************************************************************************

grid_demo: grid_search_demo.o functions.o grid_search.o global_classic.o local_classic.o
ifeq ($(__MAC__),false)
	$(CXX) $(CXXFLAGS) -o grid_demo grid_search_demo.o functions.o grid_search.o global_classic.o local_classic.o $(MPFRFLAGS)
else
	$(CXX) $(CXXFLAGS) -o grid_demo grid_search_demo.o functions.o grid_search.o global_classic.o local_classic.o $(MPFRFLAGS)
endif

#*********************************************************************************************

local_demo: local_classic_demo.o functions.o local_classic.o 
ifeq ($(__MAC__),false)
	$(CXX) $(CXXFLAGS) -o local_demo local_classic_demo.o functions.o local_classic.o $(MPFRFLAGS)
else
	$(CXX) $(CXXFLAGS) -o local_demo local_classic_demo.o functions.o local_classic.o $(MPFRFLAGS)
endif

#*********************************************************************************************

global_demo: global_classic_demo.o functions.o global_classic.o 
ifeq ($(__MAC__),false)
	$(CXX) $(CXXFLAGS) -o global_demo global_classic_demo.o functions.o global_classic.o $(MPFRFLAGS)
else
	$(CXX) $(CXXFLAGS) -o global_demo global_classic_demo.o functions.o global_classic.o $(MPFRFLAGS)
endif

#*********************************************************************************************

deflection_demo.o: deflection_demo.cpp
		$(CXX) $(CXXFLAGS) -c deflection_demo.cpp
grid_search_demo.o: grid_search_demo.cpp
		$(CXX) $(CXXFLAGS) -c grid_search_demo.cpp
local_classic_demo.o: local_classic_demo.cpp
		$(CXX) $(CXXFLAGS) -c local_classic_demo.cpp
global_classic_demo.o: global_classic_demo.cpp
		$(CXX) $(CXXFLAGS) -c global_classic_demo.cpp
functions.o: ../src/functions/functions.cpp
		$(CXX) $(CXXFLAGS) -c ../src/functions/functions.cpp
grid_search.o: ../src/grid_search/grid_search.cpp
		$(CXX) $(CXXFLAGS) -c ../src/grid_search/grid_search.cpp
global_classic.o: ../src/pso/global_classic/global_classic.cpp
		$(CXX) $(CXXFLAGS) -c ../src/pso/global_classic/global_classic.cpp
deflection.o: ../src/pso/enhanced/deflection/deflection.cpp
		$(CXX) $(CXXFLAGS) -c ../src/pso/enhanced/deflection/deflection.cpp
local_classic.o: ../src/pso/enhanced/local_classic/local_classic.cpp
		$(CXX) $(CXXFLAGS) -c ../src/pso/enhanced/local_classic/local_classic.cpp

#*********************************************************************************************

clean:
	rm -rf def_demo grid_demo local_demo global_demo *.o #grid_test