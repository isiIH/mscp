CPP=g++ -std=c++11
CPPFLAGS=-O3 -DVERBOSE
INCLUDES=-I./include/
GREEDYBW=./include/BasicCDS.cpp greedyExh.cpp
GREEDYNATURAL=./include/BasicCDS.cpp greedyExh_natural.cpp
OPTBW=./include/BasicCDS.cpp optimo.cpp
BINS=greedyExh greedyNat opt

all: clean greedyExh greedyExhNatural exhaustive

greedyExh: greedyExh.cpp
	@$(CPP) $(CPPFLAGS) $(INCLUDES) -o greedyExh $(GREEDYBW)

greedyExhNatural: greedyExh_natural.cpp
	@$(CPP) $(CPPFLAGS) $(INCLUDES) -o greedyNat $(GREEDYNATURAL)

exhaustive: optimo.cpp
	@$(CPP) $(CPPFLAGS) $(INCLUDES) -o opt $(OPTBW)

clean:
	@echo " [CLN] Removing binary files"
	@rm -f $(BINS)