# Pantranscriptome analysis C++ scripts

To compile the C++ scripts do the following (requires [CMake](https://cmake.org) 3.10 or higher): 

1. `mkdir build && cd build`
2. `cmake ..`
3. `make -j <threads>` 

Compiling the scripts should take 1-3 minutes using 4 threads (`-j`). The scripts has been successfully built on Linux (CentOS Linux 7 with GCC 8.1.0 and Ubuntu 18.04 with GCC 7.5.0) and Mac (macOS 10.14.6 with Clang 10.0.1). 

