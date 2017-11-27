APP = BPDEM2D
CC = g++
INC_PATH = ../include
LIB_PATH = ../lib

CFLAGS = -w -Wall -I$(INC_PATH)
OBJS = $(APP).o $(LIB_PATH)/*.o

all: $(APP)

$(APP).o : $(APP).cpp
	$(CC) $(CFLAGS) -c $(APP).cpp

$(APP) : $(OBJS)
	$(CC) -o $@ $(CFLAGS) $(OBJS)

clean:
	rm -f *~ \#*\#
	rm -f *.o

depend:
	makedepend -- $(CFLAGS) -- *.cpp


# DO NOT DELETE
