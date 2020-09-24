# Unmixing-SSP-Soft

1. INTRODUCTION

Contact : olivier.eches@gmail.com

This package contains programs that have been written by O. ECHES, N. DOBIGEON and J.-Y. TOURNERET, for the
implementation of the hyperspectral linear unmixing procedure for the Normal Compositional Model. For more details:

O. Eches, N. Dobigeon, C. Mailhes and J.-Y. Tourneret, "Unmixing hyperspectral images 
using a Normal Compositional Model and MCMC methods," in Proc. IEEE Workshop on Statistical 
Signal Processing (SSP), Cardiff, Wales, UK, 2009, pp.646-649.

The main library is "libUnmix" and contains the main function "unmixing.h" that calls the different sampler that 
are defined in the modelMCMC class. You can edit this file to understand the structure of the procedure based on a hierarchical model
and a Gibbs sampling strategy.
The "libUnmix" project will generate a shared library file (.so on Linux and .lib/.dll on Windows) that can be used as 3rd party library.

The program that comes with this unmixing library is the GUI "IHMUnmix" and will show you a small window where you can choose to unmix
either a real image or a synthetic pixel. Then a new window will appear with different options
for the synthetic pixel (spectral library, number of endmembers, abundance values, noise level). 
For a real image, you have to load the image data file (.tif) and the corresponding
spectral library file (.h5).

The spectral library data file for synthetic or real image must contain keys named
"Wavelength" (the wavelength number) and "spectra" (spectral values of each endmember for each wavelength number).
You can still use a hdf file viewer on the .h5 files given to see how it is constructed.

The program comes along with .tif and .h5 files for the analysis of a real image. You cand find them in the
data folder. 
Feel free to use these files as example : 
- spectraSynth.h5 for the spectral library file for a synthetic pixel
- MoffetROI.tif for the image file for unmixing of a real image
- spectraMoffetROI.h5 for the spectral library file associated with the given image file

The hyperspectral image data file must be a .tif file whose samples must be stored in a "Contiguous" configuration.
Otherwise, the program won't be able to read the file.

For the synthetic pixel unmixing option, up to 6 endmembers may be tested on a given pixel.
You can give free values to the abundance coefficients (they still must sum to one and be positive) or let
the program generate random values.
After unmixing, a window appear giving mean value of the abundances and also the histograms
of the MCMC chains of these parameters.

2. LIBRARY REQUIREMENTS

The unmixing library and the GUI programs use some libraries that must be installed before building and execution : 
- Armadillo library (http://arma.sourceforge.net/)
- Qt (http://www.qt.io)
- QCustomPlot (https://www.qcustomplot.com/index.php/introduction)
- HDF5 library (https://www.hdfgroup.org/downloads/hdf5/source-code/), it is strongly advised to download the CMake version and build the source code
- LIBTIFF with ZLIB installed (http://www.libtiff.org/)

The software has been tested under Linux Ubuntu 18.04 and under Windows 10 64bit using VS2019 (I strongly recommend to use the same software under Windows).
Please see 3.2 for how to use these libraries for building the project.

3. INSTALLATION PROCEDURE

  3.1 Procedure
  
  The "libUnmix" is using a CMakeLists.txt file for installation, simply launch CMake (feel free to use cmake-gui for better view) and select the compiler before
configuration and generation. Then, using makefile (under UNIX systems) or .sln file (under Windows with VS2019), the library can be built (do not forget the INSTALL project under VS2019 that must be built that correctly install the .lib and .dll). Use the same procedure for "IHMUnmix" but do not forget to put the generated .dll in the same folder than the executable (on Windows).

  3.2 Third party libraries

  The third party libraries should be easily available and easy to install. Some advices though :
 
 - To use the HDF5 library without errors, add HDF5_DIR as an environment variable under Linux and Windows. Add also the corresponding /bin folder to the PATH on Windows.
 - LIBTIFF and ZLIB can be installed separately. In case of difficulties, it is recommended to use the ones delivered with OTB library (https://www.orfeo-toolbox.org/download). Note that it could enter in conflict with the HDF5 library. In that case, remove all HDF5 related include files (.h) in OTB install directory, or just keep the libtiff and zlib related files. Don't forget to add to the PATH the folder where the binaries are. Note also that the OTB download files comes with Qt5 files.
 - The Armadillo library is simple to install under UNIX systems. Under Windows, you must indicate where are the include and .lib and .dll files. More precisely, the Armadillo library comes with "blas_win64_MT.lib/.dll" and "lapack_win64_MT.lib/.dll" files that must be used with any project using this library. Therefore, it must be indicated under "libUnmix/libUnmix/CMakeLists.txt" (line 33) and "IHMUnmix/CMakeLists.txt" (line 42) where to find the include and library files.
 - The QCustomPlot library can be either used "directly" (using the qcustomplot.cpp file given) or used through the shared library files that must be built. If using this last option, you will have to install QtCreator. Once built, put the generated files (debug and release) under "libUnmix/libUnmix/" directory.


