#!/bin/bash
cd doc
rm -rf doxygen-out doxybook-out website/api
doxygen Doxyfile
mkdir -p doxybook-out
doxybook2 --input doxygen-out/xml --output doxybook-out --config doxybook-config.json --templates .
# python ../../../../scripts/doc/doc_hash_replace.py ./doxybook-out/
cp -r doxybook-out website/api
cd -
