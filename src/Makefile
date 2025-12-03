CC = clang

# Compiler flags with optimizations
CFLAGS  = -Wall -Wextra -O3 -march=native -flto -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include
LDFLAGS = -L/opt/homebrew/opt/libomp/lib -flto
LDLIBS  = -lomp

# All object files
OBJECTS = acc2hit.o validate2hit.o acc2hitSerial.o

# Default target
.PHONY: all clean

all: acc2hit validate2hit acc2hitSerial

# Pattern rule for object files
%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

# Executables
acc2hit: acc2hit.o
	$(CC) -o $@ $(LDFLAGS) $^ $(LDLIBS)

validate2hit: validate2hit.o
	$(CC) -o $@ $(LDFLAGS) $^ $(LDLIBS)

acc2hitSerial: acc2hitSerial.o
	$(CC) -o $@ $(LDFLAGS) $^ $(LDLIBS)

# Dependencies (add header dependencies if any exist)
$(OBJECTS): Makefile

clean:
	rm -f *.o *~ acc2hit validate2hit acc2hitSerial

# Optional: install target
install: all
	cp acc2hit validate2hit acc2hitSerial /usr/local/bin/