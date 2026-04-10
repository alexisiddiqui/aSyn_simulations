source /homes/hussain/miniconda/bin/activate 

conda create -n gromacs python=3.11  cuda-toolkit -c conda-forge  -y
conda activate gromacs

# install dependencies
# sudo apt install cmake gcc g++ libfftw3-dev
# conda install -c conda-forge cmake fftw-devel -y


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
  -DCMAKE_CUDA_ARCHITECTURES="61;70;75;80;89" \
  -DGMX_DOUBLE=OFF \
  -DCUDA_TOOLKIT_ROOT_DIR=/homes/hussain/miniconda/envs/gromacs \
  -DCUDA_USE_STATIC_CUDA_RUNTIME=OFF \
  -DGMX_USE_HDF5=OFF \
  -DCMAKE_INSTALL_PREFIX=/ceph/opig-shared/apps/gromacs-2026

# Detect core count automatically (20 logical cores on 12700KF)
make -j$(nproc)

make check

make install