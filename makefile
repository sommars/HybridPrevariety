main:
	clear
	/bin/rm -f *.o *.out
	g++ -c polynomial_systems.cpp
	g++ -c prevariety_util.cpp
	g++ -g cone_intersection.cpp polynomial_systems.o prevariety_util.o -L/usr/include/x86_64-linux-gnu -lppl -lgmpxx -lgmp -ltbb -o prevariety.out 
	
	
# Remember to switch to g++ -O3 for both calls before benchmarking.
