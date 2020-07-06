CXX:=g++
CXXFLAGS:=-I./src -I./gfatools -O3 -std=c++14
LIBS=-lz

all: test weak_removal unbalanced_removal simple_bulge_removal bubble_removal shortcut_remover loop_killer

libgfa1.a:
	make -C gfatools
	cp gfatools/libgfa1.a .

test.o:src/test.cpp src/*.hpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

test:test.o libgfa1.a
	+$(CXX) $^ -o $@ $(LIBS)

simple_bulge_removal.o:src/simple_bulge_removal.cpp src/*.hpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

simple_bulge_removal:simple_bulge_removal.o libgfa1.a
	+$(CXX) $^ -o $@ $(LIBS)

bubble_removal.o:src/bubble_removal.cpp src/*.hpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

bubble_removal:bubble_removal.o libgfa1.a
	+$(CXX) $^ -o $@ $(LIBS)

weak_removal.o:src/weak_removal.cpp src/*.hpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

weak_removal:weak_removal.o libgfa1.a
	+$(CXX) $^ -o $@ $(LIBS)

shortcut_remover.o:src/shortcut_remover.cpp src/*.hpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

shortcut_remover:shortcut_remover.o libgfa1.a
	+$(CXX) $^ -o $@ $(LIBS)

loop_killer.o:src/loop_killer.cpp src/*.hpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

loop_killer:loop_killer.o libgfa1.a
	+$(CXX) $^ -o $@ $(LIBS)

unbalanced_removal.o:src/unbalanced_removal.cpp src/*.hpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

unbalanced_removal:unbalanced_removal.o libgfa1.a
	+$(CXX) $^ -o $@ $(LIBS)

.PHONY: clean

clean:
	rm -rf test* unbalanced_removal* loop_killer* shortcut_remover* weak_removal* bubble_removal* libgfa1.a
