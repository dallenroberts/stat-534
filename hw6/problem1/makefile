MATRICES_LSTDFLG = -lstdc++ -lm -lgsl -lgslcblas
MATRICES_INCLUDE = -I/usr/include/
MATRICES_LIB = -L/usr/lib/
MATRICES_OBJS = matrices

all:	${MATRICES_OBJS}
	rm -f *.o

main.o: main.cpp matrices.h
	gcc -g -c main.cpp -o main.o

problem1: main.o 
	gcc main.o -o matrices ${MATRICES_LIB} ${MATRICES_LSTDFLG}

clean:
	rm -f *.o
	rm -f ${MATRICES_OBJS}