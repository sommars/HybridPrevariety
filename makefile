main:
	clear
	g++ -O3 -g cone_intersection.cpp polynomial_systems.o prevariety_util.o -L/usr/include/x86_64-linux-gnu -lppl -lgmpxx -lgmp -ltbb -o prevariety.out 
	
	
# Remember to switch to g++ -O3 for both calls before benchmarking.
