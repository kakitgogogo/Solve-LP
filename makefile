CC	= g++
CFLAGS	= -std=c++11

.PHONY: test solve clean

test: test.cpp
	$(CC) $(CFLAGS) -o $@ $^

solve: solve.cpp
	$(CC) $(CFLAGS) -o $@ $^

clean: 
	rm -f test solve