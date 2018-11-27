.PHONY: all clean tests

CC=c++
CPPFLAGS=-O3 -Wall
BINARIES=bed2bedGraph bed2bedGraph-all-zero-based bed2bedGraph-all-one-based

all: $(BINARIES)

bed2bedGraph: bed2bedGraph.cpp
	$(CC) $(CPPFLAGS) -o$@ $^

bed2bedGraph-all-zero-based: bed2bedGraph-all.cpp
	$(CC) $(CPPFLAGS) -o$@ $^

bed2bedGraph-all-one-based: bed2bedGraph-all.cpp
	$(CC) -DONE_BASED $(CPPFLAGS) -o$@ $^

clean:
	rm -f $(BINARIES) 2>/dev/null && true

tests: bed2bedGraph
	cd tests; make; cd ..
