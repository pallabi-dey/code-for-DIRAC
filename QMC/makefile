CC=icx
CFLAGS=-mkl
LDFLAGS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
OBJS=multiclustU1v2.o ranlxd.o lattice.o init.o clusterops.o observable.o updatev1.o latt.o io.o

all: multiclustU1v2

multiclustU1v2: ${OBJS}
	${CC} ${CFLAGS} -o multiclustU1v2 ${OBJS} ${LDFLAGS}

clean:
	rm -rf multiclustU1v2 ${OBJS}

