CC = gcc
CFLAGS = -Wall -g -std=c99
LIBS = -lgsl -lgslcblas -lm 
INCLS = -I./ -I/usr/include/gsl


bayeslc:	bayeslc.o strparse.o matrixio.o
		${CC} ${CFLAGS} $@.o -o $@	strparse.o matrixio.o ${INCLS} ${LIBS}


.c.o: $<
		$(CC) ${INCLS} $(CFLAGS) -c $<


clean:
		\rm -f *.o *~ *%
