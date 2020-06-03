MPI_LSTDFLG = -lstdc++ -llapack -lcblas -lm
MPI_INCLUDE = -I/usr/include/
MPI_LIB = -L/usr/lib/atlas/
MPI_OBJS = parallelBayesLogistic

all:	${MPI_OBJS}
	rm -f *.o

parallelBayesLogistic.o: parallelBayesLogistic.cpp
	mpic++ -g -c parallelBayesLogistic.cpp -o parallelBayesLogistic.o ${MPI_INCLUDE}

parallelBayesLogistic: parallelBayesLogistic.o
	mpic++ parallelBayesLogistic.o -o parallelBayesLogistic ${MPI_LIB} ${MPI_LSTDFLG}

clean:
	rm -f *.o
	rm -f ${MPI_OBJS}