#!/bin/bash
cd apps/EffectiveProperties/doc/manual
rm -rf doxygen-out doxybook-out docs/manual
doxygen Doxyfile
cd doxygen-out/latex
make
cd -
cd ../../../..
