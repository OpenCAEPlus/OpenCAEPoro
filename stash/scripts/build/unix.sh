system=${1:-linux}
compiler=${2:-gnu}
build=${3:-Debug}
target=${4:-all}
if [[ ${compiler} = "intel" ]]; then
    source /opt/intel/oneapi/setvars.sh
    if [[ ${system} = "mac" ]]; then
        export CC=icc
        export CXX=icpc
    fi
fi

echo ${system}-${compiler}-${build}-${target}
cmake --preset="${system}-${compiler}-${build}" -S "."
cmake --build --preset="${system}-${compiler}-${build}" --target ${target}
