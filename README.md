Digital Signal Analysis C library (libdsac)
===========================================

Dependencies
------------

You have to install the following dependencies:

- [FFTW3](http://www.fftw.org/)
- [CMake](https://cmake.org/) >= 3.1

Quick compilation
-----------------
To compile:

```sh
git clone https://lxsd.femto-st.fr/gitlab/ahugeat/libdsac.git
cd libdsac
mkdir build
cd build
cmake ..
make
make install # may require root permissions
```

Available Options
-----------------

- DSAC_DEBUG: Enable the debug build type
- DSAC_OPENMP: Activate the OpenMP parallelization (experimental)

To use a option use this command
```sh
cmake -DDSA_DEBUG=ON|OFF ..
```
