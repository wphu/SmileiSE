/*
 * SmileiIO_Cart1D.cpp
 *
 *  Created on: 3 juil. 2013
 */

#include "SmileiIO_Cart1D.h"

#include <sstream>

#include "PicParams.h"
#include "SmileiMPI_Cart1D.h"
#include "Field1D.h"
#include "ElectroMagn.h"

using namespace std;

SmileiIO_Cart1D::SmileiIO_Cart1D( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields )
: SmileiIO( params, smpi )
{
    if(smpi->isMaster()) createPattern(params, smpi, fields);
    offset[0] = 0;
    offset[1] = 0;
    offset[2] = 0;
    offset[3] = 0;

    stride[0] = 1;
    stride[1] = 1;
    stride[2] = 1;
    stride[3] = 1;

    block[0] = 1;
    block[1] = 1;
    block[2] = 1;
    block[3] = 1;

}

SmileiIO_Cart1D::~SmileiIO_Cart1D()
{
}

//> create hdf5 data hierarchical structure: datespace, dateset and so on
void SmileiIO_Cart1D::createPattern( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields )
{
    hsize_t     dims;

    dims_global[3] = params.n_space_global[0] + 1;
    dims_global[2] = 1;
    dims_global[1] = 1;
    dims_global[0] = params.n_time / params.dump_step;

    ndims_[0] = dims_global[0];
    ndims_[1] = dims_global[1];
    ndims_[2] = dims_global[2];
    ndims_[3] = dims_global[3];


    group_id = H5Gcreate(global_file_id_, "/1d_global", H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);
    H5Gcreate(global_file_id_, "/2d_global", H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);

    /* Create a datagroup attribute. */
    dims = 4;
    //> dataspace is to descript the structure of data: the number of data dimension and the size of each dimension
    //> the first parameter rank=1: is the number of dimensions used in the dataspace
    dataspace_id = H5Screate_simple(1, &dims, NULL);
    attribute_id = H5Acreate2 (group_id, "dims_global", H5T_STD_I32BE, dataspace_id,
                                 H5P_DEFAULT, H5P_DEFAULT);
    /* Write the attribute data. */
    status = H5Awrite(attribute_id, H5T_NATIVE_INT, ndims_);
    /* Close the attribute. */
    status = H5Aclose(attribute_id);


    addField(fields->rho_global);
    addField(fields->phi_global);
    addField(fields->Ex_global);
    addField(fields->rho_global_avg);
    addField(fields->phi_global_avg);
    addField(fields->Ex_global_avg);
    for(int i = 0; i < fields->rho_s.size(); i++)
    {
        addField(fields->rho_s_global[i]);
    }

    //> if without below process, the method write() will go wrong, no ideas now!!!
    //> output initial 1d_global data===========================================
    data_ =  (double*)malloc(dims_global[3] * dims_global[2] * dims_global[1] * dims_global[0] * sizeof(double));
    for( int i = 0; i < dims_global[3] * dims_global[2] * dims_global[1] * dims_global[0]; i++)
    {
      data_[i] = 20.0;
    }
    /* Create the second dataset in group "Group_A". */
    for(int i = 0; i < dataset_id.size(); i++)
    {
        status = H5Dwrite(dataset_id[i], H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_);
    }
    free(data_);


} // END createPattern
