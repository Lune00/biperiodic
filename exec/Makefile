DEM = DEMARENA
POSTP = POSTPROCESS

CC = g++
INC_PATH = ../include
LIB_PATH = ../lib

CFLAGS = -w -Wall -I$(INC_PATH)

OBJSDEM = $(DEM).o $(LIB_PATH)/*.o
OBJSPOSTP = $(POSTP).o $(LIB_PATH)/*.o

all: $(DEM) $(POSTP)

$(DEM).o : $(DEM).cpp
	$(CC) $(CFLAGS) -c $(DEM).cpp

$(DEM) : $(OBJSDEM)
	$(CC) -o $@ $(CFLAGS) $(OBJSDEM)

$(POSTP).o : $(POSTP).cpp
	$(CC) $(CFLAGS) -c $(POSTP).cpp
$(POSTP) : $(OBJSPOSTP)
	$(CC) -o $@ $(CFLAGS) $(OBJSPOSTP)

clean:
	rm -f *~ \#*\#
	rm -f *.o

depend:
	makedepend -- $(CFLAGS) -- *.cpp


# DO NOT DELETE
