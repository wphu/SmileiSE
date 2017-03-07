/*
 * SmileiIO.cpp
 *
 *  Created on: 3 juil. 2013
 */

#include <sstream>
#include <iomanip>
#include <mpi.h>

#include "PicParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Species.h"
#include "SmileiIO.h"

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


// Create restore h5 file pattern
void SmileiIO::storeP( PicParams& params, SmileiMPI* smpi, vector<Species*>& vecSpecies, int itime )
{
    Species *s1;
    Particles *p1;
    hid_t       group_id;
    hid_t       dataspace_id;
    hid_t       dataset_id;
    hid_t       memspace_id;
    hid_t       attribute_id;
    hid_t       prop;
    herr_t      status;

    hsize_t     count[1];              /* size of subset in the file */
    hsize_t     offset[1];             /* subset offset in the file */
    hsize_t     stride[1];
    hsize_t     block[1];

    int restart;
    int partNumber;
    int nbins;
    int timestep;
    int RANK;
    double data_temp;
    hsize_t dims[1];
    hsize_t maxdims[1];
    hsize_t chunk_dims[1];

    RANK = 1;
    maxdims[0] = H5S_UNLIMITED;
    chunk_dims[0] = 2;
    restart = 0;
    timestep = 0;
    data_temp = 66.0;

    long long mpi_rk = smpi->getRank();
    string fileName = "Restore" + to_string( mpi_rk ) + "_global.h5";

    // Restore000_global.h5
    restore_file_id_  = H5Fcreate( fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    for(int iSpec = 0; iSpec < vecSpecies.size(); iSpec++)
    {
        s1 = vecSpecies[iSpec];
        p1 = &( s1->particles );
        string group_name = "/" + s1->species_param.species_type;
        group_id = H5Gcreate(restore_file_id_, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // write the variable "restart" =========================================
        dims[0]     = 1;
        count[0]    = 1;
        offset[0]   = 0;
        stride[0]   = 0;
        block[0]    = 0;
        dataspace_id = H5Screate_simple(RANK, dims, NULL);
        dataset_id = H5Dcreate2(group_id, "restart", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &restart);
        H5Sclose (dataspace_id);
        H5Dclose (dataset_id);

        // write the variable "timestep" =========================================
        dims[0]     = 1;
        count[0]    = 1;
        offset[0]   = 0;
        stride[0]   = 0;
        block[0]    = 0;
        dataspace_id = H5Screate_simple(RANK, dims, NULL);
        dataset_id = H5Dcreate2(group_id, "timestep", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &timestep);
        H5Sclose (dataspace_id);
        H5Dclose (dataset_id);

        // write the variable "PtclNumber" =================================================
        partNumber = s1->getNbrOfParticles();
        dims[0]     = 1;
        count[0]    = 1;
        offset[0]   = 0;
        stride[0]   = 0;
        block[0]    = 0;
        dataspace_id = H5Screate_simple(RANK, dims, NULL);
        dataset_id = H5Dcreate2(group_id, "partNumber", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &partNumber);
        H5Sclose (dataspace_id);
        H5Dclose (dataset_id);

        // write bmin of species ======================================================
        nbins = s1->bmin.size();
        dims[0]     = nbins;
        count[0]    = nbins;
        offset[0]   = 0;
        stride[0]   = 0;
        block[0]    = 0;
        dataspace_id = H5Screate_simple(RANK, dims, NULL);
        dataset_id = H5Dcreate2(group_id, "bmin", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(s1->bmin[0]) );
        H5Sclose (dataspace_id);
        H5Dclose (dataset_id);

        // write bmax of species ======================================================
        nbins = s1->bmax.size();
        dims[0]     = nbins;
        count[0]    = nbins;
        offset[0]   = 0;
        stride[0]   = 0;
        block[0]    = 0;
        dataspace_id = H5Screate_simple(RANK, dims, NULL);
        dataset_id = H5Dcreate2(group_id, "bmax", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(s1->bmax[0]) );
        H5Sclose (dataspace_id);
        H5Dclose (dataset_id);

        // Write momentum(0,iPart)  ======================================================
        dims[0]     = partNumber;
        count[0]    = partNumber;
        offset[0]   = 0;
        stride[0]   = 0;
        block[0]    = 0;

        dataspace_id = H5Screate_simple(RANK, dims, NULL);
        dataset_id = H5Dcreate2(group_id, "momentum0", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(p1->momentum(0,0)) );
        H5Sclose (dataspace_id);
        H5Dclose (dataset_id);


        // Write momentum(1,iPart)  ======================================================
        dims[0]     = partNumber;
        count[0]    = partNumber;
        offset[0]   = 0;
        stride[0]   = 0;
        block[0]    = 0;

        dataspace_id = H5Screate_simple(RANK, dims, NULL);
        dataset_id = H5Dcreate2(group_id, "momentum1", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(p1->momentum(1,0)) );
        H5Sclose (dataspace_id);
        H5Dclose (dataset_id);


        // Write momentum(2,iPart)  ======================================================
        dims[0]     = partNumber;
        count[0]    = partNumber;
        offset[0]   = 0;
        stride[0]   = 0;
        block[0]    = 0;

        dataspace_id = H5Screate_simple(RANK, dims, NULL);
        dataset_id = H5Dcreate2(group_id, "momentum2", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(p1->momentum(2,0)) );
        H5Sclose (dataspace_id);
        H5Dclose (dataset_id);


        // Write position(0,iPart)  ======================================================
        dims[0]     = partNumber;
        count[0]    = partNumber;
        offset[0]   = 0;
        stride[0]   = 0;
        block[0]    = 0;

        dataspace_id = H5Screate_simple(RANK, dims, NULL);
        dataset_id = H5Dcreate2(group_id, "position0", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(p1->position(0,0)) );
        H5Sclose (dataspace_id);
        H5Dclose (dataset_id);

        if(params.geometry == "2d3v")
        {
            // Write position(1,iPart)  ======================================================
            dims[0]     = partNumber;
            count[0]    = partNumber;
            offset[0]   = 0;
            stride[0]   = 0;
            block[0]    = 0;

            dataspace_id = H5Screate_simple(RANK, dims, NULL);
            dataset_id = H5Dcreate2(group_id, "position1", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status = H5Dwrite (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(p1->position(1,0)) );
            H5Sclose (dataspace_id);
            H5Dclose (dataset_id);

        }
        H5Gclose (group_id);
    }
    H5Fclose (restore_file_id_);
}

// Create restore h5 file pattern
void SmileiIO::reloadP( PicParams& params, SmileiMPI* smpi, vector<Species*>& vecSpecies, int &itime )
{
    Species *s1;
    Particles *p1;
    hid_t       group_id;
    hid_t       dataspace_id;
    hid_t       dataset_id;
    hid_t       memspace_id;
    hid_t       attribute_id;
    hid_t       prop;
    herr_t      status;

    hsize_t     count[1];              /* size of subset in the file */
    hsize_t     offset[1];             /* subset offset in the file */
    hsize_t     stride[1];
    hsize_t     block[1];

    int restart;
    int partNumber;
    int partNumber_old;
    int nbins;
    int timestep;
    int RANK;
    double data_temp;
    hsize_t dims[1];
    hsize_t dims_extented[1];
    hsize_t maxdims[1];
    hsize_t chunk_dims[1];

    RANK = 1;
    maxdims[0] = H5S_UNLIMITED;
    chunk_dims[0] = 2;
    restart = 1;
    timestep = itime;
    data_temp = 66.0;

    long long mpi_rk = smpi->getRank();
    string fileName = "Restore" + to_string( mpi_rk ) + "_global.h5";

    // Restore000_global.h5
    restore_file_id_  = H5Fopen( fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if(restore_file_id_ < 0)
    {
        return;
    }
    for(int iSpec = 0; iSpec < vecSpecies.size(); iSpec++)
    {
        s1 = vecSpecies[iSpec];
        p1 = &(s1->particles);
        string group_name = "/" + s1->species_param.species_type;
        group_id = H5Gopen(restore_file_id_, group_name.c_str(), H5P_DEFAULT);

        // Read the variable "timestep" =========================================
        dataset_id = H5Dopen(group_id, "timestep", H5P_DEFAULT);
        status = H5Dread (dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &itime);
        H5Dclose (dataset_id);

        // Read the variable "PtclNumber" =================================================
        dataset_id = H5Dopen(group_id, "partNumber", H5P_DEFAULT);
        status = H5Dread (dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &partNumber);
        H5Dclose (dataset_id);

        // resize the particles with the partNumber
        p1->initialize(partNumber, params);


        // Read momentum(0,iPart) ======================================================
        dataset_id = H5Dopen(group_id, "momentum0", H5P_DEFAULT);
        status = H5Dwrite (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(p1->momentum(0,0)) );
        H5Dclose (dataset_id);


        // Read momentum(1,iPart) ======================================================
        dataset_id = H5Dopen(group_id, "momentum1", H5P_DEFAULT);
        status = H5Dwrite (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(p1->momentum(1,0)) );
        H5Dclose (dataset_id);



        // Read momentum(2,iPart)  ======================================================
        dataset_id = H5Dopen(group_id, "momentum2", H5P_DEFAULT);
        status = H5Dwrite (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(p1->momentum(2,0)) );
        H5Dclose (dataset_id);


        // Read position(0,iPart)  ======================================================
        dataset_id = H5Dopen(group_id, "position0", H5P_DEFAULT);
        status = H5Dwrite (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(p1->position(0,0)) );
        H5Dclose (dataset_id);

        if(params.geometry == "2d3v")
        {
            // Read position(1,iPart) ======================================================
            dataset_id = H5Dopen(group_id, "position1", H5P_DEFAULT);
            status = H5Dwrite (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(p1->position(1,0)) );
            H5Dclose (dataset_id);

        }
        H5Gclose (group_id);
    }
    H5Fclose (restore_file_id_);
}




// Create restore h5 file pattern
void SmileiIO::endStoreP( PicParams& params, SmileiMPI* smpi, vector<Species*>& vecSpecies, int itime )
{
    Species *s1;
    Particles *p1;
    hid_t       group_id;
    hid_t       dataspace_id;
    hid_t       dataset_id;
    hid_t       memspace_id;
    hid_t       attribute_id;
    hid_t       prop;
    herr_t      status;

    int restart;
    restart = 0;

    long long mpi_rk = smpi->getRank();
    string fileName = "Restore" + to_string( mpi_rk ) + "_global.h5";

    // Restore000_global.h5
    restore_file_id_  = H5Fopen( fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    for(int iSpec = 0; iSpec < vecSpecies.size(); iSpec++)
    {
        s1 = vecSpecies[iSpec];
        p1 = &(s1->particles);
        string group_name = "/" + s1->species_param.species_type;
        group_id = H5Gopen(restore_file_id_, group_name.c_str(), H5P_DEFAULT);

        // write the variable "restart" =========================================
        dataset_id = H5Dopen(group_id, "restart", H5P_DEFAULT);
        status = H5Dwrite (dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &restart);
        H5Dclose (dataset_id);
        H5Gclose (group_id);
    }
    H5Fclose (restore_file_id_);
}
