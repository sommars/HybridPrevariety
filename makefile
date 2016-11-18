main:
	clear
	g++ -O3 -std=c++11 -c polynomial_systems.cpp -l:libtcmalloc_minimal.so.4.2.6
	g++ -O3 -std=c++11 -c prevariety_util.cpp -l:libtcmalloc_minimal.so.4.2.6
	g++ -O3 -std=c++11 -g -pthread cone_intersection.cpp polynomial_systems.o prevariety_util.o -l:libtcmalloc_minimal.so.4.2.6 -lppl -lgmpxx -lgmp -o prevariety.out
