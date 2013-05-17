CFLAGS = -Nmpi -fast

LFLAGS = -lpthread -lrt -lm -Nmpi

all: particle

clean:
	-$(RM) *.o particle

particle: main.o physics.o
	cc -o $@ $(CFLAGS) main.o physics.o $(LFLAGS)

touch:
	touch main.c

arc:
	tar cf - *.c *.h Makefile Makefile.sgi|gzip - > particle.tar.gz