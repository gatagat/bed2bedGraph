.PHONY: all clean tests

CC=c++
CPPFLAGS=-O3 -Wall

all: bed2bedGraph

bed2bedGraph: bed2bedGraph.cpp
	$(CC) $(CPPFLAGS) -o$@ $^

clean:
	rm -f bed2bedGraph 2>/dev/null && true

tests: bed2bedGraph
	cd tests; make; cd ..
