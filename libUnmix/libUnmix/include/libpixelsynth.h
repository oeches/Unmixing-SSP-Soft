#ifndef LIBPIXELSYNTH_H
#define LIBPIXELSYNTH_H

#include "libUnmix_global.h"
#include <iostream>
#include <armadillo>
#include "readH5File.h"

class LIBUNMIX_EXPORT libPixelSynth
{
public:
    libPixelSynth();
    ~libPixelSynth();
    int loadSpectraLib(const std::string fname);
    void generateSynthPixel(unsigned int endmNumb, unsigned int samplingStep, arma::Mat<double> abundances, double noiseVariance );
    arma::Mat<double> getSpectralValues();
    arma::Mat<double> getEmdmMat();
    arma::Mat<double> getSpectralLib();
    unsigned int getNumbBands();
    unsigned int getEndmNumb();
    unsigned int getSpectralLibrarySize();
    void setSpectralLibrarySize(unsigned int maxSize);
    void clearPixel();

private:
    unsigned int mEndmembersNumber;
    unsigned int mEndmembersNumberMax;
    arma::Mat<double> mAbMat;
    arma::Mat<double> mEndmembersMatrix;
    arma::Mat<double> mWavelength;
    arma::Mat<double> mWavelengthRed;
    arma::Mat<double> mSpectralLibrary;
    arma::Mat<double> mSpectralValues;

};

#endif // LIBPIXELSYNTH_H
