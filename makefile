CC = g++
CFLAGS = -O3 -std=c++11
ALLOC = libtcmalloc_minimal.so.4.2.6
SOPLEX = /home/jeff/Desktop/Software/soplex-2.2.1/lib/libsoplex.linux.x86_64.gnu.opt.a
PPLLINKEDFILES = -lppl -lgmpxx -lgmp -lz
OBJ = prevariety_types.o polynomial_systems.o printer.o prevariety_util.o convex_hull.o soplex_test.o process_output.o relation_tables.o cone_intersection.o

all: prevariety

prevariety: $(OBJ)
	$(CC) $(CFLAGS) -o prevariety -pthread -std=c++11 $(OBJ) $(SOPLEX) $(PPLLINKEDFILES) -l:$(ALLOC)
	
prevariety_types.o: prevariety_types.cpp prevariety_types.h
	$(CC) $(CFLAGS) -c prevariety_types.cpp

polynomial_systems.o: polynomial_systems.cpp polynomial_systems.h
	$(CC) $(CFLAGS) -c polynomial_systems.cpp

printer.o: printer.cpp printer.h
	$(CC) $(CFLAGS) -c printer.cpp

prevariety_util.o: prevariety_util.cpp prevariety_util.h
	$(CC) $(CFLAGS) -c prevariety_util.cpp

convex_hull.o: convex_hull.cpp convex_hull.h
	$(CC) $(CFLAGS) -c convex_hull.cpp

soplex_test.o: soplex_test.cpp soplex_test.h
	$(CC) $(CFLAGS) -c soplex_test.cpp

process_output.o: process_output.cpp process_output.h
	$(CC) $(CFLAGS) -c process_output.cpp

relation_tables.o: relation_tables.cpp relation_tables.h
	$(CC) $(CFLAGS) -c relation_tables.cpp
   
cone_intersection.o: cone_intersection.cpp
	$(CC) $(CFLAGS) -c cone_intersection.cpp

clean:
	/bin/rm -f *.o prevariety
