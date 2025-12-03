CC       = pgcc
CCFLAGS  = -I../common
ACCFLAGS = -acc -gpu=cuda12.3
OMPFLAGS = -fast -mp -Minfo

BIN =  acc2hit validate2hit acc2hitSerial

all: $(BIN)

acc2hit: acc2hit.c
	$(CC) $(CCFLAGS) $(ACCFLAGS) -o $@ $<

validate2hit: validate2hit.c
	$(CC) $(CCFLAGS) $(ACCFLAGS) -o $@ $<

acc2hitSerial: acc2hitSerial.c
	$(CC) $(CCFLAGS) $(ACCFLAGS) -o $@ $<

clean:
	$(RM) $(BIN)
