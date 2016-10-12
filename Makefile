CC = nvcc  -O3 -std=c++11 

CUFLAGS = -I$(CUDADIR)/include -L$(CUDADIR)/lib64  -lcuda

all: main 

main:  kernel Makefile
	$(CC) main.cpp -o main kernel.o $(CUFLAGS)  

kernel:
	$(CC) -c kernel.cu  $(CUFLAGS)
clean:
	rm  *.o main
