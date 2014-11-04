Spectral1D
==========

A simple 1D code to learn about the DFT, spectral methods, and pseudo-spectral
methods.

Available Programs
==================

- transform: simply transforms a signal to/from spectral space.
- differentiate: differentiates a signal using spectral means.
- aliasing: examines aliasing errors arising from nonlinear product terms.
- advection: solves the linear advection equation with a spectral approach for
  the spatial derivative and a TVD RK3 for time integration.
- advectionFD: solves the linear advection equation with a 2nd order FD
  approach for the spatial derivative and a TVD RK3 for time integration.

Required Software and Environment Variables
===========================================

- You must have FFTW installed (v3.3.3 was used when writing this package).
- You must have the environment variable FFTW_ROOT set to the root directory
  of the FFTW installation. For example, if your FFTW libraries and header
  files are installed in /usr/lib and /usr/include, set FFTW_ROOT to /usr.
  In bash, placing "export FFTW_ROOT=/usr" in your .bashrc will work.
- You must have cmake version 2.8.9 or higher installed.

Installation
============

1. Change to the main source directory:
    $ cd /path/to/Spectral1D
2. Make a build directory and go into it:
    $ mkdir build
    $ cd build
3. Enter the cmake curses GUI and press c to configure cmake:
    $ ccmake ..
    $ [in cmake curses gui] c
4. Toggle which of the programs you would like to build:
    $ [in cmake curses gui] enter
5. Enter the installation prefix (note these can run from build just fine):
    $ [in cmake curses gui] CMAKE_INSTALL_PREFIX=/path/to/installation
6. Configure your changes in cmake:
    $ [in cmake curses gui] c
7. Generate the Makefile and exit cmake curses gui:
    $ [in cmake curses gui] g
8. Make the programs.
    $ make
9. Install the programs.
    $ make install
