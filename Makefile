CPP=g++ -std=c++11
CPPFLAGS=-O3 -DVERBOSE
INCLUDES=-I./include/
GREEDY=./include/BasicCDS.cpp greedy.cpp
GREEDYBW=./include/BasicCDS.cpp greedyBW.cpp
OPTBW=./include/BasicCDS.cpp optimo.cpp
BINS=greedy greedyBW opt

all: clean greedy greedyBW exhaustive

greedy: greedy.cpp
	@$(CPP) $(CPPFLAGS) $(INCLUDES) -o greedy $(GREEDY)

greedyBW: greedyBW.cpp
	@$(CPP) $(CPPFLAGS) $(INCLUDES) -o greedyBW $(GREEDYBW)

exhaustive: optimo.cpp
	@$(CPP) $(CPPFLAGS) $(INCLUDES) -o opt $(OPTBW)

clean:
	@echo " [CLN] Removing binary files"
	@rm -f $(BINS)