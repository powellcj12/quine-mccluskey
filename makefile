CC = g++
CFLAGS = -c -Wall -Wextra
all: QM

QM: QM.o
	$(CC) QM.o -o QM

QM.o: QM.cpp
	$(CC) $(CFLAGS) QM.cpp

clean:
	rm -rf *o QM