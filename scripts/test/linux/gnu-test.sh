cmake -G "Ninja" -DCMAKE_BUILD_TYPE="Release" -DCMAKE_CXX_COMPILER="g++" -DCMAKE_C_COMPILER="gcc" -DCMAKE_Fortran_COMPILER="gfortran" -DFASP_DIR="/opt/faspsolver" -DBUILD_TEST=ON -S "." -B "./build/x64-gnu-Release"
cmake --build "./build/x64-gnu-Release" --target all
ctest --test-dir "./build/x64-gnu-Release"