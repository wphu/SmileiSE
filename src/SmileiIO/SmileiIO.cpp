/*
 * SmileiIO.cpp
 *
 *  Created on: 3 juil. 2013
 */

#include "SmileiIO.h"

#include <sstream>
#include <iomanip>

#include <mpi.h>

#include "PicParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Species.h"

using namespace std;

// static varable must be defined and initialized here

SmileiIO::SmileiIO( PicParams& params, SmileiMPI* smpi )
{
    // Fields_global.h5
    if(smpi->isMaster()) global_file_id_  = H5Fcreate( "Fields_global.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    ndims_t = 0;

}

SmileiIO::~SmileiIO()
{
    // Management of global IO file
    //H5Fclose( global_file_id_ );
}




void SmileiIO::addField(Field* field)
{
    string fieldName = field->name;
    const char* name = fieldName.c_str();
    dataset_name.push_back(name);

    /* Create the data space for the dataset. */
    dataspace_id = H5Screate_simple(4, dims_global, NULL);
    int dataset_size = dataset_id.size();
    hid_t id = H5Dcreate2(group_id, name, H5T_NATIVE_DOUBLE, dataspace_id,
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    dataset_id.push_back(id);

    dataset_field.push_back(field);
}

//! write potential, rho and so on into hdf5 file every some timesteps
void SmileiIO::write(  ElectroMagn* fields, SmileiMPI* smpi )
{
    offset[0] = ndims_t;


    count[0]  = 1;
    count[1]  = dims_global[1];
    count[2]  = dims_global[2];
    count[3]  = dims_global[3];


    //>write the ===global potential=== in time ndims_t to the file
    for(int i = 0; i < dataset_id.size(); i++)
    {
        memspace_id = H5Screate_simple (4, count, NULL);
        dataspace_id = H5Dget_space (dataset_id[i]);
        status = H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset,
                                          stride, count, block);
        status = H5Dwrite (dataset_id[i], H5T_NATIVE_DOUBLE, memspace_id,
                             dataspace_id, H5P_DEFAULT, dataset_field[i]->data_);

        status = H5Sclose (memspace_id);
        status = H5Sclose (dataspace_id);
    }

    ndims_t++;

} // END write
