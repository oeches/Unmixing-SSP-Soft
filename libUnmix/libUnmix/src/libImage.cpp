#include <iostream>
#include "include/libImage.h"
#include "tiffio.h"
#include <armadillo>


LibImage::LibImage()
{
    mDimensionX = 1;
    mDimensionY = 1;
    mSpectralBands = 1;
    mEndmembersNumber = 1;
    mWavelength = arma::zeros(mSpectralBands,1);
    mEndmembersMatrix = arma::zeros(mSpectralBands,mEndmembersNumber);
    mPixelValues = arma::zeros(mDimensionX*mDimensionY,mSpectralBands);
    mSpectralValues = arma::zeros(mSpectralBands,1);
}

LibImage::~LibImage()
{

}

int LibImage::loadImage(const std::string fnameImage, const std::string fnameData)
{
    //reading image
    TIFF* tif = TIFFOpen(fnameImage.c_str(), "r");
    int retReadHD = -1;
    arma::Mat<uint8_t> MatMask;
    try
    {
        if (tif)
        {
            // reading spectra
            retReadHD = readHDF5ToMat<double>(fnameData, "lib dataset", "spectra", mEndmembersMatrix);

            // reading wavelength
            retReadHD = readHDF5ToMat<double>(fnameData, "lib dataset", "Wavelength", mWavelength);

            // reading mask
            retReadHD = readHDF5ToMat<uint8_t>(fnameData, "lib dataset", "mask", MatMask);
        }
        else
        {
            throw std::string("Error : image file don't exists !");
        }
    }
    catch (std::string &str)
    {
        std::cerr << str << std::endl;
        return EXIT_FAILURE;
    }
    // putting on arma matrix
    if (retReadHD == 0)
    {
        // extracting pixel
        short config = 0;
        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &mDimensionX);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &mDimensionY);
        TIFFGetField(tif, TIFFTAG_PLANARCONFIG, &config);
        uint32 imagelength;
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &imagelength);

        //get number of bands
        arma::uvec vecMask = arma::find(MatMask);
        mSpectralBands = vecMask.n_rows;
        uint32 numBandsRaw = MatMask.n_cols;
        mEndmembersNumber = mEndmembersMatrix.n_cols;

        // tif buffer
        tsize_t sls = TIFFScanlineSize(tif);
        tdata_t buf = _TIFFmalloc(sls);

        //intermediate arma matrix
        arma::mat newbuf;
        try
        {
            if (config == PLANARCONFIG_CONTIG)
            {
                newbuf.set_size(mDimensionX * numBandsRaw, mDimensionY);

                for (uint32 row = 0; row < imagelength; ++row)
                {
                    TIFFReadScanline(tif, buf, row, 0);
                    memcpy(newbuf.colptr(row), buf, sls);
                }
                _TIFFfree(buf);
                buf = nullptr;
                TIFFClose(tif);

                mPixelValues.set_size(mDimensionX * mDimensionY, mSpectralBands);

                for (uint32 k = 0; k < mDimensionY; ++k)
                {
                    arma::mat interBuf;
                    interBuf.set_size(mDimensionX, numBandsRaw);
                    interBuf = arma::trans(arma::reshape(newbuf.col(k), numBandsRaw, mDimensionX));
                    uint32 bandNum = 0;
                    for (uint32 i = 0; i < numBandsRaw; ++i)
                    {
                        if (MatMask(0, i) == 1)
                        {
                            mPixelValues(arma::span(k * mDimensionX, mDimensionX * (k + 1) - 1), bandNum) = interBuf.col(i);
                            ++bandNum;
                        }
                    }
                }
            }
            else
            {
                // NOT SUPPORTED
                throw std::string("Error : current config of file not supported !");
            }
        }
        catch (std::string& strConf)
        {
            std::cerr << strConf << std::endl;
            return EXIT_FAILURE;
        }
    }

    return retReadHD;
}

void LibImage::getPixel(const unsigned int index)
{
    mSpectralValues.set_size(1,mSpectralBands);
    mSpectralValues = mPixelValues.row(index);
    mSpectralValues.reshape(mSpectralBands,1);
}

arma::Mat<double> LibImage::getSpectralValues()
{
    return mSpectralValues;
}

arma::Mat<double> LibImage::getEmdmMat()
{
    return mEndmembersMatrix;
}

unsigned int LibImage::getEndmNumb()
{
    return mEndmembersNumber;
}

unsigned int LibImage::getNPixels()
{
    return mDimensionX*mDimensionY;
}

unsigned int LibImage::getDimX()
{
    return mDimensionX;
}

unsigned int LibImage::getDimY()
{
    return mDimensionY;
}
