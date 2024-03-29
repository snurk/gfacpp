CXX:=g++
CXXFLAGS:=-I./src -I./gfatools -I./gfakluge -O3 -std=c++14 -Wall
LIBS:=-lz
ODIR:=build
DEPS:=src/*.hpp
#SRCS=$(wildcard src/*.cpp)
#EXECS=$(patsubst src/%.cpp,$(ODIR)/%,$(SRCS))
EXECS:=test neighborhood unambig_extension weak_removal unbalanced_removal simple_bulge_removal bubble_removal shortcut_remover loop_killer nongenomic_link_removal tip_clipper low_cov_remover isolated_remover

all: $(patsubst %,$(ODIR)/%,$(EXECS))

#all: $EXECS

$(ODIR)/libgfa1.a:gfatools/*.c gfatools/*.h
	make -C gfatools
	cp gfatools/libgfa1.a $@

$(ODIR)/wrapper.o:src/wrapper.cpp $(DEPS)
	$(CXX) -c $(CXXFLAGS) $< -o $@

$(ODIR)/%.o:src/%.cpp $(DEPS)
	$(CXX) -c $(CXXFLAGS) $< -o $@

$(ODIR)/%:$(ODIR)/%.o $(ODIR)/wrapper.o $(ODIR)/libgfa1.a
	$(CXX) $^ -o $@ $(LIBS)

.PRECIOUS: $(ODIR)/%.o

.PHONY: clean
clean:
	rm -rf $(ODIR)/*
	rm -f gfatools/libgfa1.a
	make -C gfatools clean
