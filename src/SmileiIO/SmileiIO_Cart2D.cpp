/*
 * SmileiIO_Cart2D.cpp
 *
 *  Created on: 3 juil. 2013
 */

#include "SmileiIO_Cart2D.h"

#include <sstream>

#include "PicParams.h"
#include "SmileiMPI_Cart2D.h"
#include "Field2D.h"
#include "ElectroMagn.h"

using namespace std;

SmileiIO_Cart2D::SmileiIO_Cart2D( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies )
: SmileiIO( params, smpi )
{
    if(smpi->isMaster())
    {
        createPattern(params, smpi, fields);
        status = H5Fclose(global_file_id_);
    }
    fieldsGroup.offset[0] = 0;
    fieldsGroup.offset[1] = 0;
    fieldsGroup.offset[2] = 0;
    fieldsGroup.offset[3] = 0;

    fieldsGroup.stride[0] = 1;
    fieldsGroup.stride[1] = 1;
    fieldsGroup.stride[2] = 1;
    fieldsGroup.stride[3] = 1;

    fieldsGroup.block[0] = 1;
    fieldsGroup.block[1] = 1;
    fieldsGroup.block[2] = 1;
    fieldsGroup.block[3] = 1;

}

SmileiIO_Cart2D::~SmileiIO_Cart2D()
{
}

//> create hdf5 file, datespace, dateset and so on
void SmileiIO_Cart2D::createPattern( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields )
{
    hsize_t     dims;

    fieldsGroup.dims_global[3] = params.n_space_global[1] + 1;
    fieldsGroup.dims_global[2] = params.n_space_global[0] + 1;
    fieldsGroup.dims_global[1] = 1;
    fieldsGroup.dims_global[0] = params.n_time / params.dump_step;

    fieldsGroup.ndims_[0] = fieldsGroup.dims_global[0];
    fieldsGroup.ndims_[1] = fieldsGroup.dims_global[1];
    fieldsGroup.ndims_[2] = fieldsGroup.dims_global[2];
    fieldsGroup.ndims_[3] = fieldsGroup.dims_global[3];


    data_ =  (double*)malloc(fieldsGroup.dims_global[3] * fieldsGroup.dims_global[2] * fieldsGroup.dims_global[1] * fieldsGroup.dims_global[0] * sizeof(double));
    for( int i = 0; i < fieldsGroup.dims_global[3] * fieldsGroup.dims_global[2] * fieldsGroup.dims_global[1] * fieldsGroup.dims_global[0]; i++)
    {
      data_[i] = 20.0;
    }

    fieldsGroup.group_id = H5Gcreate(global_file_id_, "/Fields", H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);

    /* Create a datagroup attribute. */
    dims = 4;
    fieldsGroup.dataspace_id = H5Screate_simple(1, &dims, NULL);
    fieldsGroup.attribute_id = H5Acreate2 (fieldsGroup.group_id, "dims_global", H5T_STD_I32BE, fieldsGroup.dataspace_id,
                                 H5P_DEFAULT, H5P_DEFAULT);
    /* Write the attribute data. */
    fieldsGroup.status = H5Awrite(fieldsGroup.attribute_id, H5T_NATIVE_INT, fieldsGroup.ndims_);
    /* Close the attribute. */
    fieldsGroup.status = H5Aclose(fieldsGroup.attribute_id);


    addField(fields->rho_global);
    addField(fields->phi_global);
    addField(fields->Ex_global);
    addField(fields->rho_global_avg);
    addField(fields->phi_global_avg);
    addField(fields->Ex_global_avg);

    //> if without below process, the method write() will go wrong, no ideas now!!!
    //> output initial 1d_global data===========================================
    data_ =  (double*)malloc(fieldsGroup.dims_global[3] * fieldsGroup.dims_global[2] * fieldsGroup.dims_global[1] * fieldsGroup.dims_global[0] * sizeof(double));
    for( int i = 0; i < fieldsGroup.dims_global[3] * fieldsGroup.dims_global[2] * fieldsGroup.dims_global[1] * fieldsGroup.dims_global[0]; i++)
    {
      data_[i] = 20.0;
    }
    /* Create the second dataset in group "Group_A". */
    for(int i = 0; i < fieldsGroup.dataset_id.size(); i++)
    {
        fieldsGroup.status = H5Dwrite(fieldsGroup.dataset_id[i], H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_);
    }
    free(data_);


} // END createPattern
