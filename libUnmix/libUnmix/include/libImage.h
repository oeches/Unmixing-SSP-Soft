#ifndef LIBIMAGE_H
#define LIBIMAGE_H

#include "libUnmix_global.h"
#include <iostream>
#include <armadillo>
#include "readH5File.h"

class LIBUNMIX_EXPORT LibImage
{
public:
    LibImage();
    ~LibImage();
    int loadImage(const std::string fnameImage, const std::string fnameData);
    void getPixel(const unsigned int index);
    arma::Mat<double> getSpectralValues();
    arma::Mat<double> getEmdmMat();
    unsigned int getEndmNumb();
    unsigned int getNPixels();
    unsigned int getDimX();
    unsigned int getDimY();

private:
    unsigned int mDimensionX;
    unsigned int mDimensionY;
    unsigned int mSpectralBands;
    unsigned int mEndmembersNumber;
    arma::Mat<double> mWavelength;
    arma::Mat<double> mEndmembersMatrix;
    arma::Mat<double> mPixelValues;
    arma::Mat<double> mSpectralValues;
};



#endif // LIBIMAGE_H
