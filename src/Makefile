CC = gcc
CFLAGS = -Wall -O2
LDLIBS =
OBJECTS = acc2hit.o validate2hit.o

all: acc2hit validate2hit

acc2hit.o: acc2hit.c

validate2hit.o: validate2hit.c

acc2hit: acc2hit.o
	$(CC) -o acc2hit $(CFLAGS) $(LDLIBS) acc2hit.o

validate2hit: validate2hit.o
	$(CC) -o validate2hit $(CFLAGS) $(LDLIBS) validate2hit.o

clean:
	rm -f *.o *~ acc2hit validate2hit