RND_LSTDFLG = -lstdc++ -lm -lgsl -lgslcblas
RND_INCLUDE = -I/usr/include/
RND_LIB = -L/usr/lib/
RND_OBJS = rnd

all:	${RND_OBJS}
	rm -f *.o

matrices.o: matrices.cpp matrices.h
	gcc -g -c matrices.cpp -o matrices.o ${RND_INCLUDE}

main.o: main.cpp matrices.h
	gcc -g -c main.cpp -o main.o

rnd: main.o
	gcc main.o matrices.o -o rnd ${RND_LIB} ${RND_LSTDFLG}

clean:
	rm -f *.o
	rm -f ${RND_OBJS}