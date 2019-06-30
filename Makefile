CC=g++
CFLAGS=-c -Wall -I ~/boost_1_55_0 -O3
#CFLAGS=-c -Wall -I /opt/boost -O3
LDFLAGS=
SOURCES=main.cc model.cc solution.cc loopupdate.cc workflow.cc measure.cc
OBJECTS=$(SOURCES:.cc=.o)
EXECUTABLE=bilayer

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cc.o:
	$(CC) $(CFLAGS) $< -o $@