APP = BPDEM2D
DEM = DEMARENA
CC = g++
INC_PATH = ../include
LIB_PATH = ../lib

CFLAGS = -w -Wall -I$(INC_PATH)
OBJS = $(APP).o $(LIB_PATH)/*.o

OBJSDEM = $(DEM).o $(LIB_PATH)/*.o

all: $(APP) $(DEM)

$(APP).o : $(APP).cpp
	$(CC) $(CFLAGS) -c $(APP).cpp

$(APP) : $(OBJS)
	$(CC) -o $@ $(CFLAGS) $(OBJS)

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
