#include "include/libpixelsynth.h"
#include <iostream>
#include <armadillo>

libPixelSynth::libPixelSynth()
{
    mEndmembersNumber = 1;
    mEndmembersNumberMax = 6;
    mAbMat = arma::zeros(mEndmembersNumber,1);
    mWavelength = arma::zeros(1,1);
    mWavelengthRed = arma::zeros(1,1);
    mEndmembersMatrix = arma::zeros(1,mEndmembersNumber);
    mSpectralLibrary = arma::zeros(1, mEndmembersNumberMax);
    mSpectralValues = arma::zeros(1,1);
}

libPixelSynth::~libPixelSynth()
{
    
}

int libPixelSynth::loadSpectraLib(const std::string fname)
{
    // putting on arma matrix
    int retReadHD;
    // reading spectra
    retReadHD= readHDF5ToMat<double>(fname, "lib dataset", "spectra", mSpectralLibrary);

    if (retReadHD == 0)
    {
        mEndmembersNumberMax = mSpectralLibrary.n_cols;
    }

    // reading wavelength
    retReadHD= readHDF5ToMat<double>(fname, "lib dataset", "Wavelength", mWavelength);

    return retReadHD;
}

void libPixelSynth::generateSynthPixel(unsigned int endmNumb, unsigned int samplingStep, arma::Mat<double> abundances, double noiseVariance )
{
    unsigned int numberofBandsOld, numberofBandsOldRed;
    arma::vec meanVec;

    mEndmembersNumber = endmNumb;
    mAbMat.set_size(mEndmembersNumber,1);
    if (arma::size(abundances) == arma::size(mAbMat))
    {
        mAbMat = abundances;
    }
    else
    {
        throw "Error : wrong number of endmembers";
    }

    arma::vec indxEndmChosenVec = arma::linspace(0, endmNumb-1, endmNumb);
    arma::uvec indxEndmChosenUVec = arma::conv_to<arma::uvec>::from(indxEndmChosenVec);
    numberofBandsOld = mWavelength.n_elem;
    numberofBandsOldRed = (unsigned int) numberofBandsOld/samplingStep;
    arma::vec indxBandsVec = arma::linspace(0, numberofBandsOld-1, numberofBandsOldRed);
    arma::uvec indxBandsUVec = arma::conv_to<arma::uvec>::from(indxBandsVec);

    mEndmembersMatrix.set_size(numberofBandsOldRed, endmNumb);
    mEndmembersMatrix = mSpectralLibrary.submat(indxBandsUVec, indxEndmChosenUVec);
    mWavelengthRed.set_size(numberofBandsOldRed,1);
    mWavelengthRed.rows(0,numberofBandsOldRed-1) = mWavelength.rows(indxBandsUVec);

    arma::Mat<double> sigMat = noiseVariance* arma::eye(numberofBandsOldRed,numberofBandsOldRed);
    arma::Mat<double> randMatrixEnd(numberofBandsOldRed, endmNumb);

    for (unsigned int i=0; i < endmNumb; ++i)
    {
        meanVec = mEndmembersMatrix.col(i);
        randMatrixEnd.col(i) = mvnrnd(meanVec, sigMat);
    }
    mSpectralValues.set_size(numberofBandsOldRed, 1);
    mSpectralValues = randMatrixEnd*mAbMat;
}

arma::Mat<double> libPixelSynth::getSpectralValues()
{
    return mSpectralValues;
}

arma::Mat<double> libPixelSynth::getEmdmMat()
{
    return mEndmembersMatrix;
}

arma::Mat<double> libPixelSynth::getSpectralLib()
{
    return mSpectralLibrary;
}

unsigned int libPixelSynth::getNumbBands()
{
    return mWavelength.n_rows;
}

unsigned int libPixelSynth::getEndmNumb()
{
    return mEndmembersNumber;
}

unsigned int libPixelSynth::getSpectralLibrarySize()
{
    return mEndmembersNumberMax;
}

void libPixelSynth::setSpectralLibrarySize(unsigned int maxSize)
{
    mEndmembersNumberMax = maxSize;
    mSpectralLibrary = mSpectralLibrary.cols(0, maxSize-1);
}

void libPixelSynth::clearPixel()
{
    mSpectralValues = arma::zeros(1,1);
    mEndmembersMatrix = arma::zeros(1,mEndmembersNumber);
    mAbMat = arma::zeros(mEndmembersNumber,1);
    mWavelengthRed = arma::zeros(1,1);
}
