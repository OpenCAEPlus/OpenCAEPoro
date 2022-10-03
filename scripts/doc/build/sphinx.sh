#!/bin/bash
cd doc
make latexpdf
cp _build/latex/opencaeporo.pdf website/manual/
cd -
