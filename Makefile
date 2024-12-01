CPP=g++ -std=c++11
CPPFLAGS=-O3 -DVERBOSE -fopenmp
INCLUDES=-I./include/
GRASP=./include/BasicCDS.cpp ./src/graspSC.cpp
GRASPCPU=./include/BasicCDS.cpp ./src/graspSC_CPU.cpp
# OPTCU=./include/BasicCDS.cpp optimo_cuda.cu
BINS=grasp graspCPU

all: clean grasp grasp_cpu

grasp: src/graspSC.cpp
	@$(CPP) $(CPPFLAGS) $(INCLUDES) -o grasp $(GRASP)

grasp_cpu: src/graspSC_CPU.cpp
	@$(CPP) $(CPPFLAGS) $(INCLUDES) -o graspCPU $(GRASPCPU)

# optimo_cuda: optimo_cuda.cu
# 	nvcc $(CUDAFLAGS) $(INCLUDES) -o opt_cu $(OPTCU)

%.o: %.cu
	nvcc $(CUDAFLAGS) -c $< -o $@

clean:
	@echo " [CLN] Removing binary files"
	@rm -f $(BINS)