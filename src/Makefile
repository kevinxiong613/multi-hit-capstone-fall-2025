CC = /opt/homebrew/opt/llvm/bin/clang

CFLAGS  = -Wall -O2 -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include
LDFLAGS = -L/opt/homebrew/opt/libomp/lib
LDLIBS  = -lomp

OBJECTS = acc2hit.o validate2hit.o

all: acc2hit validate2hit

acc2hit.o: acc2hit.c

validate2hit.o: validate2hit.c

acc2hit: acc2hit.o
	$(CC) -o acc2hit $(CFLAGS) $(LDFLAGS) acc2hit.o $(LDLIBS)

validate2hit: validate2hit.o
	$(CC) -o validate2hit $(CFLAGS) $(LDFLAGS) validate2hit.o $(LDLIBS)

clean:
	rm -f *.o *~ acc2hit validate2hit
