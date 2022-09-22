cmake -G "Ninja" -DCMAKE_BUILD_TYPE="$1" -DCMAKE_CXX_COMPILER="g++" -DCMAKE_C_COMPILER="gcc" -DCMAKE_Fortran_COMPILER="gfortran" -DFASP_DIR="/opt/faspsolver" -S "." -B "./build/x64-gnu-$1"
cmake --build "./build/x64-gnu-$1" --target OpenCAEPoro
