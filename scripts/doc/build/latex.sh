#!/bin/bash
cd doc
rm -rf doxygen-out
doxygen Doxyfile
cd doxygen-out/latex
make
cd -
cd ..
