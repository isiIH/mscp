CPP=g++ -std=c++11
CPPFLAGS=-O3 -DVERBOSE -fopenmp
INCLUDES=-I./include/
OPTBW=./include/BasicCDS.cpp optimo.cpp
BINS=opt

all: clean optimo

optimo: optimo.cpp
	@$(CPP) $(CPPFLAGS) $(INCLUDES) -o opt $(OPTBW)

clean:
	@echo " [CLN] Removing binary files"
	@rm -f $(BINS)