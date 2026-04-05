curl -O https://ftp.gromacs.org/gromacs/gromacs-2026.1.tar.gz
tar xfz gromacs-2026.1.tar.gz
cd gromacs-2026.1
mkdir build && cd build

export PATH="/opt/homebrew/opt/llvm/bin:$PATH"
export CC=/opt/homebrew/opt/llvm/bin/clang
export CXX=/opt/homebrew/opt/llvm/bin/clang++
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=/opt/homebrew/opt/llvm/bin/clang \
  -DCMAKE_CXX_COMPILER=/opt/homebrew/opt/llvm/bin/clang++ \
  -DOpenMP_C_FLAGS="-fopenmp" \
  -DOpenMP_CXX_FLAGS="-fopenmp" \
  -DOpenMP_C_LIB_NAMES="omp" \
  -DOpenMP_CXX_LIB_NAMES="omp" \
  -DOpenMP_omp_LIBRARY=/opt/homebrew/opt/libomp/lib/libomp.dylib \
  -DGMX_SIMD=ARM_NEON_ASIMD \
  -DGMX_FFT_LIBRARY=fftw3 \
  -DGMX_BUILD_OWN_FFTW=ON \
  -DGMX_OPENMP=ON \
  -DGMX_THREAD_MPI=ON \
  -DGMX_GPU=OpenCL \
  -DGMX_DOUBLE=OFF \
  -DCMAKE_INSTALL_PREFIX=$HOME/gromacs-2026


# Use all cores (-j$(sysctl -n hw.logicalcpu) detects core count automatically)
make -j$(sysctl -n hw.logicalcpu)

# Run regression tests (optional but recommended)
make check

# Install
make install