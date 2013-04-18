
CFLAGS = -Nmpi

LFLAGS= -lpthread -lrt -lm

CC= icc

all: blurc thresc

clean:
	-$(RM) *.o blur thresc

blurc: ppmio.o gaussw.o blurfilter.o blurmain.o
	cc -o $@ ppmio.o gaussw.o blurfilter.o blurmain.o $(LFLAGS)

thresc: thresmain.o ppmio.o thresfilter.o
	icc -Nmpi -o $@ thresmain.o ppmio.o thresfilter.o $(LFLAGS)

arc:
	tar cf - *.c *.h *.f90 Makefile Makefile.sgi|gzip - > filters.tar.gz
