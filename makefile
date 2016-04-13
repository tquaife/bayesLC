CC = gcc
HOME=/home/tqu
CFLAGS = -Wall -g -std=c99
#CFLAGS = -Wall -O3
#CFLAGS = -Wall -std=c99
LIBS =  -L${HOME}/lib -lmatrixio -lstrparse  -lgsl -lgslcblas -lm 

STATIC = -static -L${HOME}/lib ${HOME}/lib/libstrparse.a ${HOME}/lib/libmatrixio.a /usr/lib/libgsl.a /usr/lib/libgslcblas.a /usr/lib/libm.a

INCLS = -I./ -I${HOME}/include -I/usr/include/gsl
BINDIR = ${HOME}/bin/linux

VER=bayeslc

bayeslc:	bayeslc.o
		${CC} ${CFLAGS}  $@.o -o $@	${INCLS} ${LIBS}
		#mv bayeslc ${BINDIR}

static:		${VER}.o
		${CC} -o ${VER} ${VER}.o ${CFLAGS} ${INCLS} ${STATIC}  	


.c.o: $<
		$(CC) ${INCLS} $(CFLAGS) -c $<


clean:
		\rm -f *.o *~ *%
