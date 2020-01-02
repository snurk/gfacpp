CXX:=g++
CXXFLAGS:=-I./src -I./gfatools -O3 -std=c++14
LIBS=-lz

all: test weak_removal

libgfa1.a:
	make -C gfatools
	cp gfatools/libgfa1.a .

test.o:src/test.cpp src/wrapper.hpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

test:test.o libgfa1.a
	+$(CXX) $^ -o $@ $(LIBS)

weak_removal.o:src/weak_removal.cpp src/wrapper.hpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

weak_removal:weak_removal.o libgfa1.a
	+$(CXX) $^ -o $@ $(LIBS)

.PHONY: clean

clean:
	rm -rf test* libgfa1.a
