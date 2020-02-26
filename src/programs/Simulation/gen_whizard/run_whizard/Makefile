CC=g++
CFLAGS= -c -O3 -g -Wall `root-config --cflags`
LDFLAGS= -g -O3 -D_REENTRANT `root-config --glibs` -lTreePlayer
SOURCES=run_whizard.cc MyReadConfig.cc
OBJECTS=$(SOURCES:.cc=.o)
EXECUTABLE=run_whizard

all: $(SOURCES) $(EXECUTABLE)
	

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cc.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f ./*~ ./*.o ./run_whizard
