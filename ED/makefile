CC=icx
CXX=icpx
EIGEN_DIR=/usr/include/eigen3
LIBS=-lm -L/opt/intel/oneapi/mkl/latest/lib/intel64 -qmkl -L/usr/include/boost/   
OBJS=common.o construct.o Hamiltonian.o observable.o io.o 2dU1v2.o
CFLAGS= -std=c99 -Wall -qopenmp
CXXFLAGS= -std=c++11 -Wall -mcmodel=large -DMKL_ILP64 -I$(EIGEN_DIR) -qopenmp

all: 2dU1v2

2dU1v2: ${OBJS}
	${CXX} -o 2dU1v2 ${OBJS} ${LIBS}
clean:
	rm -rf 2dU1v2 ${OBJS} *.optrpt 
	
	
	
	
	

