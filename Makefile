CPP=g++ -std=c++11
CPPFLAGS=-O3 -DVERBOSE
INCLUDES=-I./include/
GREEDYBW=./include/BasicCDS.cpp greedyExh.cpp
OPTBW=./include/BasicCDS.cpp optimo.cpp
BINS=greedyExh opt

all: clean greedyExh exhaustive

greedyExh: greedyExh.cpp
	@$(CPP) $(CPPFLAGS) $(INCLUDES) -o greedyExh $(GREEDYBW)

exhaustive: optimo.cpp
	@$(CPP) $(CPPFLAGS) $(INCLUDES) -o opt $(OPTBW)

clean:
	@echo " [CLN] Removing binary files"
	@rm -f $(BINS)