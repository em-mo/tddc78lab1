CFLAGS = -Nmpi -O2
CXXFLAGS = -Nmpi -O2

LFLAGS = -lpthread -lrt -lm -Nmpi

CC = icc
CXX = icpc

all: particle

clean:
	-$(RM) *.o particle

particle: main.o physics.o
	icc -o $@ $(CFLAGS) main.o physics.o $(LFLAGS)

touch:
	touch main.cpp

arc:
	tar cf - *.cpp *.h Makefile Makefile.sgi|gzip - > particle.tar.gz