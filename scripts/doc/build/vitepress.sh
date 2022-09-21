#!/bin/bash
pnpm --filter opencaeporo-doc doc:build
cd doc/website/.vitepress
zip -r dist.zip dist
cd -
