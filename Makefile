CPP = g++ -std=c++17
CPPFLAGS = -O3 -g -DVERBOSE -fopenmp
INCLUDES = -I./include/
LDFLAGS = -ltbb

SRC_DIR = ./include

SRC_MAIN = $(SRC_DIR)/BasicCDS.cpp $(SRC_DIR)/Set.cpp $(SRC_DIR)/SCP.cpp $(SRC_DIR)/Group.cpp \
      $(SRC_DIR)/Grasp.cpp $(SRC_DIR)/SetCover.cpp $(SRC_DIR)/RowCovering.cpp \
      ./main.cpp

BINS=grasp

all: clean $(BINS)

grasp: main.cpp
	@$(CPP) $(CPPFLAGS) $(INCLUDES) -o grasp $(SRC_MAIN) $(LDFLAGS)

clean:
	@echo " [CLN] Removing binary files"
	@rm -f $(BINS)

# for ((i=0; i<=5; i++)); do ./grasp scp $((2**i)); done;