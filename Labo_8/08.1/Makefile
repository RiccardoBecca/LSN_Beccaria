CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o lib.o funzionebase.o integral.o
	$(CC) random.o lib.o main.o funzionebase.o integral.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
lib.o : lib.cpp lib.h
	$(CC) -c lib.cpp -o lib.o $(CFLAGS)
funzionebase.o : funzionebase.cpp funzionebase.h
	$(CC) -c funzionebase.cpp -o funzionebase.o $(CFLAGS)
integral.o : integral.cpp integral.h
	$(CC) -c integral.cpp -o integral.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
