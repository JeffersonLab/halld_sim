CC=g++
CFLAGS= -c -O3 -g -Wall `root-config --cflags`
LDFLAGS= -g -O3 -D_REENTRANT `root-config --glibs` -lTreePlayer
SOURCES=lhe_to_root.cc 
OBJECTS=$(SOURCES:.cc=.o)
EXECUTABLE=lhe_to_root

all: $(SOURCES) $(EXECUTABLE)
	

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cc.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f ./*~ ./*.o ./lhe_to_root
