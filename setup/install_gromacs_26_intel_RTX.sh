sudo apt install cmake gcc g++ libfftw3-dev

# Download & extract (same as before)
curl -O https://ftp.gromacs.org/gromacs/gromacs-2026.1.tar.gz
tar xfz gromacs-2026.1.tar.gz
cd gromacs-2026.1
mkdir build && cd build



# No special compiler overrides needed on Linux — GCC works great
# Make sure you have: sudo apt install cmake gcc g++ libfftw3-dev cuda-toolkit-XX
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DGMX_SIMD=AVX2_256 \
  -DGMX_FFT_LIBRARY=fftw3 \
  -DGMX_BUILD_OWN_FFTW=ON \
  -DGMX_OPENMP=ON \
  -DGMX_THREAD_MPI=ON \
  -DGMX_GPU=CUDA \
  -DGMX_CUDA_TARGET_SM=86 \
  -DGMX_DOUBLE=OFF \
  -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda \
  -DCUDA_USE_STATIC_CUDA_RUNTIME=OFF \
  -DCMAKE_INSTALL_PREFIX=$HOME/gromacs-2026

# Detect core count automatically (20 logical cores on 12700KF)
make -j$(nproc)

make check

make install