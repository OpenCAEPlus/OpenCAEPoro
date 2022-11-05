#!/bin/bash
# Regression test for OpenCAEPoro project

clear
date '+Regression test at: %A, %B %d, %Y %H:%M:%S'

RED='\033[0;31m'   # red color
GREEN='\033[0;32m' # green color
NC='\033[0m'       # no color

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

PNUM_ALL=8
if (($# > 0)) && (($1<=PNUM_ALL)); then
    PNUM=$1
else
    PNUM=$PNUM_ALL
fi
printf "Runing ${RED}${PNUM}/${PNUM_ALL}${NC} regression tests ...\n"

OCP_LOG="Log.out"
OCP_SUM="Summary.out"
OCP_REV="FastReview.out"

for ((i=0; i<PNUM; i++)); do
    if [ ! -d "$i" ]; then
        mkdir $i
    fi

    printf "Problem $i: ${GREEN}${DIRS[i]}/${FILES[i]}${NC}\n"
    ../testOpenCAEPoro ../examples/${DIRS[i]}/${FILES[i]} > ../examples/${DIRS[i]}/log.out
    mv    ../examples/${DIRS[i]}/$OCP_LOG  ./$i/New_IMPEC_$OCP_LOG
    mv    ../examples/${DIRS[i]}/$OCP_SUM  ./$i/New_IMPEC_$OCP_SUM
    mv    ../examples/${DIRS[i]}/$OCP_REV  ./$i/New_IMPEC_$OCP_REV
    sdiff  ./$i/Old_IMPEC_$OCP_LOG ./$i/New_IMPEC_$OCP_LOG | grep "|" 2>/dev/null
    sdiff  ./$i/Old_IMPEC_$OCP_SUM ./$i/New_IMPEC_$OCP_SUM | grep "|" 2>/dev/null
    sdiff  ./$i/Old_IMPEC_$OCP_REV ./$i/New_IMPEC_$OCP_REV | grep "|" 2>/dev/null
done

# Overwrite results
OVERWRITE_OLD="N"
printf "${RED}Attention${NC}: Overwrite the existing results? [yes/N]: "
read OVERWRITE_OLD

if [ "$OVERWRITE_OLD" = "yes" ]; then
    for ((i=0; i<PNUM; i++)); do
        mv $i/New_FIM_$OCP_LOG   $i/Old_FIM_$OCP_LOG   2>/dev/null
        mv $i/New_FIM_$OCP_SUM   $i/Old_FIM_$OCP_SUM   2>/dev/null
        mv $i/New_FIM_$OCP_REV   $i/Old_FIM_$OCP_REV   2>/dev/null
        mv $i/New_IMPEC_$OCP_LOG $i/Old_IMPEC_$OCP_LOG 2>/dev/null
        mv $i/New_IMPEC_$OCP_SUM $i/Old_IMPEC_$OCP_SUM 2>/dev/null
        mv $i/New_IMPEC_$OCP_REV $i/Old_IMPEC_$OCP_REV 2>/dev/null
    done
    echo "Overwrite the old results."
else
    echo "Keep both old and new results."
fi
