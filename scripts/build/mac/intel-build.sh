source /opt/intel/oneapi/setvars.sh
CC=icc CXX=icpc cmake -G "Ninja" -DFASP_DIR="/opt/faspsolver" -DCMAKE_BUILD_TYPE="$1" -DCMAKE_CXX_COMPILER="icpc" -DCMAKE_C_COMPILER="icc" -S "." -B "./build/x64-intel-$1"
cmake --build "./build/x64-intel-$1" --target all