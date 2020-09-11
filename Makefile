CXX?=g++
CXXFLAGS:=-I./src -I./gfatools -I./gfakluge -O3 -std=c++14
LIBS:=-lz
ODIR:=build
DEPS:=src/*.hpp
#SRCS=$(wildcard src/*.cpp)
#EXECS=$(patsubst src/%.cpp,$(ODIR)/%,$(SRCS))
EXECS:=test weak_removal unbalanced_removal simple_bulge_removal bubble_removal shortcut_remover loop_killer

LEGACY_DEPS:=src/legacy/*.hpp gfakluge/*.hpp
LEGACY_EXECS:=neighborhood cnt_aware_tip_clipper tip_clipper isolated_remover low_cov_remover path_length #bubble_finder

all: $(patsubst %,$(ODIR)/%,$(EXECS)) $(patsubst %,$(ODIR)/%,$(LEGACY_EXECS))

#all: $EXECS

$(ODIR)/libgfa1.a:gfatools/*.c gfatools/*.h
	make -C gfatools
	cp gfatools/libgfa1.a $@

$(ODIR)/%.o:src/%.cpp $(DEPS)
	$(CXX) -c $(CXXFLAGS) $< -o $@

$(ODIR)/%.o:src/legacy/%.cpp $(DEPS) $(LEGACY_DEPS)
	$(CXX) -c $(CXXFLAGS) $< -o $@

$(ODIR)/%:$(ODIR)/%.o $(ODIR)/libgfa1.a
	$(CXX) $^ -o $@ $(LIBS)

.PRECIOUS: $(ODIR)/%.o

.PHONY: clean
clean:
	rm -rf $(ODIR)/*
