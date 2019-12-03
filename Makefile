CXX:=g++
CXXFLAGS:=-I./src -I./gfatools -O3 -std=c++14
LIBS=-lz

all: test

libgfa1.a:
	make -C gfatools
	cp gfatools/libgfa1.a .

test.o:src/test.cpp src/wrapper.hpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

test:test.o libgfa1.a
	+$(CXX) $^ -o $@ $(LIBS)

.PHONY: clean

clean:
	rm -rf test* libgfa1.a
