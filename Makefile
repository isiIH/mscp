CPP=g++ -std=c++11
CPPFLAGS=-O3 -DVERBOSE -fopenmp
INCLUDES=-I./include/
OPTBW=./include/BasicCDS.cpp optimo.cpp
OPTCU=./include/BasicCDS.cpp optimo_cuda.cu
BINS=opt opt_cu

CUDAFLAGS  := -arch=sm_50

all: clean optimo optimo_cuda

optimo: optimo.cpp
	@$(CPP) $(CPPFLAGS) $(INCLUDES) -o opt $(OPTBW)

optimo_cuda: optimo_cuda.cu
	nvcc $(CUDAFLAGS) $(INCLUDES) -o opt_cu $(OPTCU)

%.o: %.cu
	nvcc $(CUDAFLAGS) -c $< -o $@

clean:
	@echo " [CLN] Removing binary files"
	@rm -f $(BINS)