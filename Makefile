CPP=g++ -std=c++11
CPPFLAGS=-O3 -DVERBOSE -fopenmp
INCLUDES=-I./include/
GRASP=./include/BasicCDS.cpp graspSC.cpp
# OPTCU=./include/BasicCDS.cpp optimo_cuda.cu
BINS=grasp

all: clean grasp

grasp: graspSC.cpp
	@$(CPP) $(CPPFLAGS) $(INCLUDES) -o grasp $(GRASP)

# optimo_cuda: optimo_cuda.cu
# 	nvcc $(CUDAFLAGS) $(INCLUDES) -o opt_cu $(OPTCU)

%.o: %.cu
	nvcc $(CUDAFLAGS) -c $< -o $@

clean:
	@echo " [CLN] Removing binary files"
	@rm -f $(BINS)