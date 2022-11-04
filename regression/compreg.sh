#!/bin/bash
# Regression test for OpenCAEPoro project

clear
date '+Regression test at: %A, %B %d, %Y %H:%M:%S'

RED='\033[0;31m'   # red color
GREEN='\033[0;32m' # green color
NC='\033[0m'       # no color

OCP_LOG="log.out"
OCP_SUM="SUMMARY.out"
OCP_REV="FastReview.out"

PNUM=8
DIRS=("spe1a"
      "spe1a"
      "spe1b"
      "spe1b" 
      "spe9"
      "cornerpoint"
      "spe3"
      "spe5")
FILES=("spe1a.data" 
       "spe1a.data method=FIM dtInit=1 dtMax=100 dtMin=0.1"
       "spe1b.data"
       "spe1b.data method=FIM dtInit=1 dtMax=100 dtMin=0.1" 
       "spe9_FIM.data"
       "CP.data    method=FIM dtInit=1 dtMax=100 dtMin=0.1"
       "spe3.data  method=FIM dtInit=1 dtMax=100 dtMin=0.1"
       "spe5.data  method=FIM dtInit=1 dtMax=50  dtMin=0.1")

for ((i=0; i<PNUM; i++)); do
    if [ ! -d "$i" ]; then
        mkdir $i
    fi

    printf "Solve ${GREEN}${DIRS[i]}/${FILES[i]}${NC} ...\n"
    ../testOpenCAEPoro ../examples/${DIRS[i]}/${FILES[i]} > ../examples/${DIRS[i]}/log.out
    mv    ../examples/${DIRS[i]}/$OCP_LOG  ./$i/New_IMPEC_$OCP_LOG
    mv    ../examples/${DIRS[i]}/$OCP_SUM  ./$i/New_IMPEC_$OCP_SUM
    mv    ../examples/${DIRS[i]}/$OCP_REV  ./$i/New_IMPEC_$OCP_REV
    sdiff  ./$i/Old_IMPEC_$OCP_LOG ./$i/New_IMPEC_$OCP_LOG | grep "|"
    sdiff  ./$i/Old_IMPEC_$OCP_SUM ./$i/New_IMPEC_$OCP_SUM | grep "|"
    sdiff  ./$i/Old_IMPEC_$OCP_REV ./$i/New_IMPEC_$OCP_REV | grep "|"

done

# Overwrite results
OVERWRITE_OLD="N"
echo "Attention: Overwrite the existing results? [yes/N]: "
read OVERWRITE_OLD

if [ "$OVERWRITE_OLD" = "yes" ]; then
    for ((i=0; i<PNUM; i++)); do
        mv $i/New_FIM_$OCP_LOG   $i/Old_FIM_$OCP_LOG
        mv $i/New_FIM_$OCP_SUM   $i/Old_FIM_$OCP_SUM
        mv $i/New_FIM_$OCP_REV   $i/Old_FIM_$OCP_REV
        mv $i/New_IMPEC_$OCP_LOG $i/Old_IMPEC_$OCP_LOG
        mv $i/New_IMPEC_$OCP_SUM $i/Old_IMPEC_$OCP_SUM
        mv $i/New_IMPEC_$OCP_REV $i/Old_IMPEC_$OCP_REV
    done
    echo "Overwrite the existing results."
else
    echo "Keep old and new results."
fi
