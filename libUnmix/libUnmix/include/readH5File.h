#ifndef READH5FILE_H
#define READH5FILE_H

#include "H5Cpp.h"
#include <armadillo>
#include <vector>

using namespace H5;

template<typename T>
int readHDF5ToMat(H5std_string fname, H5std_string groupname, H5std_string dataset, arma::Mat<T> &dataMatOut);

template <typename T>
class DataTypeFor
{
public:
};


template <>
class DataTypeFor<uint8_t>
{
public:
    const H5::DataType value = H5::PredType::NATIVE_UINT8;
};

template <>
class DataTypeFor<double>
{
public:
    const H5::DataType value = H5::PredType::NATIVE_DOUBLE;
};

template<typename T>
int readHDF5ToMat(H5std_string fname, H5std_string groupname, H5std_string dataset, arma::Mat<T> &dataMatOut)
{
    try
    {
        H5::Exception::dontPrint();

        // Open an existing file, group and dataset.
        H5File file(fname, H5F_ACC_RDONLY );
        Group group = file.openGroup(groupname);
        DataSet datasetDat = group.openDataSet(dataset);

        H5::DataType dt;
        DataTypeFor<T> valueDat;
        dt = valueDat.value;
        DataSpace dataspace = datasetDat.getSpace();
        hsize_t dims_out[2];
        int ndims = dataspace.getSimpleExtentDims(dims_out);
        if (ndims == 1)
        {
            dims_out[1] = 1;
        }
        //std::vector<T> data_out;
        T* data_out = new T[dims_out[0]*dims_out[1]];
        datasetDat.read(data_out, dt);
        dataMatOut.resize(dims_out[0], dims_out[1]);
        dataMatOut.fill(arma::fill::zeros);

        for (hsize_t i = 0; i < dims_out[0]; ++i)
        {
            for (hsize_t j = 0; j < dims_out[1]; ++j)
            {
                dataMatOut(i,j) = data_out[j+i*dims_out[1]];
            }
        }
        if (data_out !=0)
        {
            delete[] data_out;
            data_out = 0;
        }
        return 0;

    }

    // catch failure caused by the H5File operations
    catch(H5::FileIException error)
    {
        error.printErrorStack();
        return -1;
    }

    // catch failure caused by the DataSet operations
    catch(H5::DataSetIException error)
    {
        error.printErrorStack();
        return -1;
    }

}

#endif // READH5FILE_H
