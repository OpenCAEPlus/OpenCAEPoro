#!/bin/bash
cd doc
rm -rf doxybook-out website/api
mkdir -p doxybook-out
doxybook2 --input doxygen-out/xml --output doxybook-out --config doxybook-config.json --templates .
cp -r doxybook-out website/api
cd -
