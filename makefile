CC = g++
CFLAGS = -O3 -std=c++11
ALLOC = libtcmalloc_minimal.so.4.2.6
PPLLINKEDFILES = -lppl -lgmpxx -lgmp
OBJ = polynomial_systems.o prevariety_util.o cone_intersection.o

all: prevariety

prevariety: $(OBJ)
	$(CC) $(CFLAGS) -o prevariety -pthread -std=c++11 $(OBJ) $(PPLLINKEDFILES) -l:$(ALLOC)
	
polynomial_systems.o: polynomial_systems.cpp polynomial_systems.h
	$(CC) $(CFLAGS) -c polynomial_systems.cpp

prevariety_util.o: prevariety_util.cpp prevariety_util.h
	$(CC) $(CFLAGS) -c prevariety_util.cpp

cone_intersection.o: cone_intersection.cpp
	$(CC) $(CFLAGS) -c cone_intersection.cpp

clean:
	/bin/rm -f *.o prevariety
