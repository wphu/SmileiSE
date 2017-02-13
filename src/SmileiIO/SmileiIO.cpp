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
