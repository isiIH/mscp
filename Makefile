CPP = g++ -std=c++17
CPPFLAGS = -O3 -DVERBOSE -fopenmp
INCLUDES = -I./include/

SRC_DIR = ./include

SRC_MAIN = $(SRC_DIR)/BasicCDS.cpp $(SRC_DIR)/Set.cpp $(SRC_DIR)/SCP.cpp $(SRC_DIR)/UnionFind.cpp \
      $(SRC_DIR)/Grasp.cpp $(SRC_DIR)/SetCover.cpp $(SRC_DIR)/RowCovering.cpp $(SRC_DIR)/Group.cpp \
      ./main.cpp

# GRASPCPU=$(SRC_DIR)/BasicCDS.cpp ./graspSC_CPU.cpp

BINS=grasp # graspCPU

all: clean $(BINS)

grasp: main.cpp
	@$(CPP) $(CPPFLAGS) $(INCLUDES) -o grasp $(SRC_MAIN)

# graspCPU: graspSC_CPU.cpp
# 	@$(CPP) $(CPPFLAGS) $(INCLUDES) -o graspCPU $(GRASPCPU)

clean:
	@echo " [CLN] Removing binary files"
	@rm -f $(BINS)