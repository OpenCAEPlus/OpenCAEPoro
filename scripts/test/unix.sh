system=${1:-linux}
compiler=${2:-gnu}
build=${3:-Debug}
source ./scripts/build/unix.sh ${system} ${compiler} ${build}
ctest --preset="${system}-${compiler}-${build}"
