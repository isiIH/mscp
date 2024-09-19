CPP=g++ -std=c++11
CPPFLAGS=-O3 -DVERBOSE -fopenmp
INCLUDES=-I./include/
APRBW=./include/BasicCDS.cpp aproxSC.cpp
# OPTCU=./include/BasicCDS.cpp optimo_cuda.cu
BINS=apr

all: clean aprox

aprox: aproxSC.cpp
	@$(CPP) $(CPPFLAGS) $(INCLUDES) -o apr $(APRBW)

# optimo_cuda: optimo_cuda.cu
# 	nvcc $(CUDAFLAGS) $(INCLUDES) -o opt_cu $(OPTCU)

%.o: %.cu
	nvcc $(CUDAFLAGS) -c $< -o $@

clean:
	@echo " [CLN] Removing binary files"
	@rm -f $(BINS)