# Overview

Library for common operations on GFA graphs and some processing algorithms.
*wrapper.hpp* provides a C++ interface from gfatools (https://github.com/lh3/gfatools/) implementation of sequence graphs and some additional functionality.
Legacy algorithms (in src/legacy folder) were implemented based on a gfakluge (https://github.com/edawson/gfakluge) graph implementation (included in gfakluge folder) and will be reimplemented via the gfatools wrapper.

# Get and compile

To clone:
```
git clone --recurse-submodules git@github.com:snurk/gfacpp.git
```

To compile:
```
cd gfacpp ; make
```

Each binary in build folder corresponds to one processing algorithm.
Run binary required parameters.

# Description of individual procedures
TBD

# TODO

* provide documentation
* add description of individual procedures
* implement graph compaction
* implement unambigous_traversal.cpp
* migrate to clipp.hpp
