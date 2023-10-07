CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o lib.o TSP_GA.o
	$(CC) TSP_GA.o random.o lib.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
lib.o : ../nifty_lib/lib.cpp ../nifty_lib/lib.h
	$(CC) -c ../nifty_lib/lib.cpp -o lib.o $(CFLAGS)
TSP_GA.o : TSP_GA.cpp TSP_GA.h
	$(CC) -c TSP_GA.cpp -o TSP_GA.o $(CFLAGS)
	
clean :
	rm *.o main.exe seed.out
