DEM = DEMARENA
CC = g++
INC_PATH = ../include
LIB_PATH = ../lib

CFLAGS = -w -Wall -I$(INC_PATH)

OBJSDEM = $(DEM).o $(LIB_PATH)/*.o

all: $(DEM)

$(DEM).o : $(DEM).cpp
	$(CC) $(CFLAGS) -c $(DEM).cpp

$(DEM) : $(OBJSDEM)
	$(CC) -o $@ $(CFLAGS) $(OBJSDEM)

clean:
	rm -f *~ \#*\#
	rm -f *.o

depend:
	makedepend -- $(CFLAGS) -- *.cpp


# DO NOT DELETE
