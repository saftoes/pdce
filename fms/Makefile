INCS=-I/usr/local/include -I/usr/local/include/flint
LIBS=-L/usr/local/lib -larb -lflint -lmpirxx -lmpir -lmpfr -lgmp -lm -lpthread 

CC=gcc
CFLAGS=-ansi -pedantic -Wall -O2 -funroll-loops -g -mpopcnt 

all:
	$(CC) $(INCS) $(CFLAGS) -o main.x main.cpp $(LIBS)