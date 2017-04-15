CC	= g++
CFLAGS	= -std=c++11

.PHONY: test

test: test.cpp
	$(CC) $(CFLAGS) -o $@ $^
