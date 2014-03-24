CC=g++
CFLAGS=-c -Wall -m64
LDFLAGS=
SOURCES=main.cpp field.cpp geometry.cpp grid.cpp gridarea.cpp human.cpp maplocation.cpp random.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=MobilitySimulator

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
clean:
	rm -rf *o MobilitySimulator
