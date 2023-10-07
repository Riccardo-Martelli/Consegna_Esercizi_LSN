CC = g++
CFLAGS = -Wall -O3 --std=c++11

NVE_NVT.exe : MD_MC.o random.o lib.o
	$(CC) random.o MD_MC.o lib.o -o NVE_NVT.exe
MD_MC.o : MD_MC.cpp
	$(CC) -c MD_MC.cpp -o MD_MC.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
lib.o : ../../../nifty_lib/lib.cpp ../../../nifty_lib/lib.h
	$(CC) -c ../../../nifty_lib/lib.cpp -o lib.o $(CFLAGS)
clean :
	rm *.o NVE_NVT.exe
