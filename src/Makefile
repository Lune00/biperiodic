INC_PATH = ../include
LIB_PATH = ../lib
CC = g++
CFLAGS = -w -Wall -I$(INC_PATH)
SRC = Algo.cpp Particle.cpp Cell.cpp Config.cpp Tenseur.cpp Vecteur.cpp Sample.cpp Interactions.cpp Contact.cpp globalfunctions.cpp

OBJS = $(subst .cpp,.o,$(SRC))
.SUFFIXES: .cpp
.cpp.o:
	$(CC) -c $(CFLAGS) $<
	cp -f $@ $(LIB_PATH)

all: $(OBJS)

clean:
	rm -f *~ \#*\#
	rm -f *.o

depend:
	makedepend -- $(CFLAGS) -- *.cpp

cleandist: clean
	makedepend clean
