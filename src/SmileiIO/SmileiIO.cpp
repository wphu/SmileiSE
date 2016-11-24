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
    fieldsGroup.dataset_name.push_back(name);

    /* Create the data space for the dataset. */
    fieldsGroup.dataspace_id = H5Screate_simple(4, fieldsGroup.dims_global, NULL);
    int dataset_size = fieldsGroup.dataset_id.size();
    hid_t id = H5Dcreate2(fieldsGroup.group_id, name, H5T_NATIVE_DOUBLE, fieldsGroup.dataspace_id,
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    fieldsGroup.dataset_id.push_back(id);

    fieldsGroup.dataset_data.push_back(field->data_);
}

//! write potential, rho and so on into hdf5 file every some timesteps
void SmileiIO::write( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies)
{

    calVDF( params, smpi, fields, vecSpecies);

    if(smpi->isMaster()) {
        fieldsGroup.offset[0] = ndims_t;
        ptclsGroup.offset[0] = ndims_t;

        fieldsGroup.count[0]  = 1;
        fieldsGroup.count[1]  = fieldsGroup.dims_global[1];
        fieldsGroup.count[2]  = fieldsGroup.dims_global[2];
        fieldsGroup.count[3]  = fieldsGroup.dims_global[3];

        ptclsGroup.count[0]  = 1;
        ptclsGroup.count[1]  = ptclsGroup.dims_global[1];
        ptclsGroup.count[2]  = ptclsGroup.dims_global[2];
        ptclsGroup.count[3]  = ptclsGroup.dims_global[3];

        //>write fields
        for(int i = 0; i < fieldsGroup.dataset_id.size(); i++)
        {
            fieldsGroup.memspace_id = H5Screate_simple (4, fieldsGroup.count, NULL);
            fieldsGroup.dataspace_id = H5Dget_space (fieldsGroup.dataset_id[i]);
            fieldsGroup.status = H5Sselect_hyperslab (fieldsGroup.dataspace_id, H5S_SELECT_SET, fieldsGroup.offset,
                                              fieldsGroup.stride, fieldsGroup.count, fieldsGroup.block);
            fieldsGroup.status = H5Dwrite (fieldsGroup.dataset_id[i], H5T_NATIVE_DOUBLE, fieldsGroup.memspace_id,
                                 fieldsGroup.dataspace_id, H5P_DEFAULT, fieldsGroup.dataset_data[i]);

            fieldsGroup.status = H5Sclose (fieldsGroup.memspace_id);
            fieldsGroup.status = H5Sclose (fieldsGroup.dataspace_id);
        }

        // write particle velocity distribution function
        for(int i = 0; i < ptclsGroup.dataset_id.size(); i++)
        //for(int i = 0; i < 0; i++)
        {
            ptclsGroup.memspace_id = H5Screate_simple (4, ptclsGroup.count, NULL);
            ptclsGroup.dataspace_id = H5Dget_space (ptclsGroup.dataset_id[i]);
            ptclsGroup.status = H5Sselect_hyperslab (ptclsGroup.dataspace_id, H5S_SELECT_SET, ptclsGroup.offset,
                                              ptclsGroup.stride, ptclsGroup.count, ptclsGroup.block);
            ptclsGroup.status = H5Dwrite (ptclsGroup.dataset_id[i], H5T_NATIVE_DOUBLE, ptclsGroup.memspace_id,
                                 ptclsGroup.dataspace_id, H5P_DEFAULT, ptclsGroup.dataset_data[i]);

            ptclsGroup.status = H5Sclose (fieldsGroup.memspace_id);
            ptclsGroup.status = H5Sclose (fieldsGroup.dataspace_id);
        }


        ndims_t++;
    }


} // END write
