CXX:=g++
CXXFLAGS:=-I./src -I./gfatools -O3 -std=c++14

all: test

libgfa1.a:
	make -C gfatools
	cp gfatools/libgfa1.a .

test.o:src/test.cpp src/wrapper.hpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

test:test.o libgfa1.a
	+$(CXX) -o $@ $^

.PHONY: clean

clean:
	rm -rf test* libgfa1.a
