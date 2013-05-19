CFLAGS = -Nmpi -fast

LFLAGS = -lpthread -lrt -lm -Nmpi

CC = icc

all: particle

clean:
	-$(RM) *.o particle

particle: main.o physics.o
	icc -o $@ $(CFLAGS) main.o physics.o $(LFLAGS)

touch:
	touch main.cpp

arc:
	tar cf - *.cpp *.h Makefile Makefile.sgi|gzip - > particle.tar.gz