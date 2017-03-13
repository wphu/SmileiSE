/*
 * SmileiIO.cpp
 *
 *  Created on: 3 juil. 2013
 */
#include <sys/io.h>
#include <sys/dir.h>

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
    is_restart = 0;
    stepStart  = 0;
    global_file_name_ = "data_global.h5";
    restore_file_dir_ = "restore";

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
    fieldsGroup.dataset_data.push_back(field->data_);
}


void SmileiIO::createFieldsGroup( ElectroMagn* fields )
{

    // ============ Create fieldsGroup ============================================
    //addField(fields->rho_global);
    //addField(fields->phi_global);
    //addField(fields->Ex_global);
    addField(fields->rho_global_avg);
    addField(fields->phi_global_avg);
    addField(fields->Ex_global_avg);
    for(int i = 0; i < fields->rho_s.size(); i++)
    {
        //addField(fields->rho_s_global[i]);
        addField(fields->rho_s_global_avg[i]);
        addField(fields->Vx_s_global_avg[i]);
        addField(fields->Vy_s_global_avg[i]);
        addField(fields->Vz_s_global_avg[i]);
        addField(fields->Vp_s_global_avg[i]);
        addField(fields->T_s_global_avg[i]);

    }

    if(!is_restart)
    {
        fieldsGroup.group_id = H5Gcreate(global_file_id_, "/Fields", H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);
        // ============Create hdf5 struct for fieldsGroup===================================
        for(int i = 0; i < fieldsGroup.dataset_name.size(); i++)
        {
            // Create hdf5 datasets under "/Fields" group
            fieldsGroup.dataspace_id = H5Screate_simple(4, fieldsGroup.dims_global, NULL);
            hid_t id = H5Dcreate2(fieldsGroup.group_id, fieldsGroup.dataset_name[i], H5T_NATIVE_DOUBLE, fieldsGroup.dataspace_id,
                                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            fieldsGroup.dataset_id.push_back(id);
            H5Dclose( fieldsGroup.dataset_id[i] );
            H5Sclose( fieldsGroup.dataspace_id );
        }

        // ==============Create and write the attribute data ================================
        //> dataspace is to descript the structure of data: the number of data dimension and the size of each dimension
        //> the first parameter rank=1: is the number of dimensions used in the dataspace
        fieldsGroup.dataspace_id = H5Screate_simple(1, &(fieldsGroup.aDims), NULL);
        fieldsGroup.attribute_id = H5Acreate2 (fieldsGroup.group_id, "dims_global", H5T_STD_I32BE, fieldsGroup.dataspace_id,
                                     H5P_DEFAULT, H5P_DEFAULT);
        fieldsGroup.status = H5Awrite(fieldsGroup.attribute_id, H5T_NATIVE_INT, fieldsGroup.ndims_);
        H5Aclose( fieldsGroup.attribute_id );
        H5Sclose( fieldsGroup.dataspace_id );


        H5Gclose( fieldsGroup.group_id );
    }
    else
    {
        fieldsGroup.dataset_id.resize( fieldsGroup.dataset_name.size() );
    }
}

void SmileiIO::createPartsGroup( vector<Species*>& vecSpecies )
{
    Species *s;
    Particles *p;
    // ==============Create "VDF" group ================================
    for(int isp=0; isp<vx_VDF.size(); isp++)
    {
        s = vecSpecies[isp];

        string fieldName = "VDF_" + s->species_param.species_type;
        ptclsGroup.dataset_stringName.push_back(fieldName);
        const char* name = ptclsGroup.dataset_stringName.back().c_str();
        ptclsGroup.dataset_name.push_back(name);
        ptclsGroup.dataset_data.push_back(vx_VDF_global[isp]->data_);

        // for VDF_tot_ the dimension is same with VDF_, but only [nx=0] is meaningful
        fieldName = "VDF_tot_" + s->species_param.species_type;
        ptclsGroup.dataset_stringName.push_back(fieldName);
        name = ptclsGroup.dataset_stringName.back().c_str();
        ptclsGroup.dataset_name.push_back(name);
        ptclsGroup.dataset_data.push_back(vx_VDF_tot_global[isp]->data_);

    }


    if(!is_restart)
    {
        ptclsGroup.group_id = H5Gcreate(global_file_id_, "/VDF", H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);
        // ==============Create hdf5 struct for "VDF" group ================================
        for(int i = 0; i < ptclsGroup.dataset_name.size(); i++)
        {
            ptclsGroup.dataspace_id = H5Screate_simple(4, ptclsGroup.dims_global, NULL);
            hid_t id = H5Dcreate2(ptclsGroup.group_id, ptclsGroup.dataset_name[i], H5T_NATIVE_DOUBLE, ptclsGroup.dataspace_id,
                                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            ptclsGroup.dataset_id.push_back(id);
            H5Dclose( ptclsGroup.dataset_id[i] );
            H5Sclose( ptclsGroup.dataspace_id );
        }

        // ==============Create and write the attribute data ================================
        //> dataspace is to descript the structure of data: the number of data dimension and the size of each dimension
        //> the first parameter rank=1: is the number of dimensions used in the dataspace
        ptclsGroup.dataspace_id = H5Screate_simple(1, &(ptclsGroup.aDims), NULL);
        ptclsGroup.attribute_id = H5Acreate2 (ptclsGroup.group_id, "dims_global", H5T_STD_I32BE, ptclsGroup.dataspace_id,
                                     H5P_DEFAULT, H5P_DEFAULT);
        ptclsGroup.status = H5Awrite(ptclsGroup.attribute_id, H5T_NATIVE_INT, ptclsGroup.ndims_);
        // close attribute, dataset, dataspace, group, and so on
        H5Aclose( ptclsGroup.attribute_id );
        H5Sclose( ptclsGroup.dataspace_id );

        H5Gclose( ptclsGroup.group_id );
    }
    else
    {
        ptclsGroup.dataset_id.resize( ptclsGroup.dataset_name.size() );
    }
}

// Create restore h5 file pattern
void SmileiIO::storeP( PicParams& params, SmileiMPI* smpi, vector<Species*>& vecSpecies, int itime )
{
    Species     *s1;
    Particles   *p1;
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
    restart = 1;
    timestep = itime;
    data_temp = 66.0;
    nSpecies = vecSpecies.size();

    long long mpi_rk = smpi->getRank();
    string fileName = restore_file_dir_ + "/Restore" + to_string( mpi_rk ) + "_global.h5";

    // Restore000_global.h5
    restore_file_id_  = H5Fcreate( fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // ========================== Write info about restart ============================
    string group_name = "/info" ;
    group_id = H5Gcreate(restore_file_id_, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // write the variable "restart"
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

    // write the variable "timestep"
    dims[0]     = 1;
    count[0]    = 1;
    offset[0]   = 0;
    stride[0]   = 0;
    block[0]    = 0;
    dataspace_id = H5Screate_simple(RANK, dims, NULL);
    dataset_id = H5Dcreate2(group_id, "timestep", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite (dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &timestep);
    H5Sclose( dataspace_id );
    H5Dclose( dataset_id );

    // write the variable "nSpecies"
    dims[0]     = 1;
    count[0]    = 1;
    offset[0]   = 0;
    stride[0]   = 0;
    block[0]    = 0;
    dataspace_id = H5Screate_simple(RANK, dims, NULL);
    dataset_id = H5Dcreate2(group_id, "nSpecies", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite (dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nSpecies);
    H5Sclose( dataspace_id );
    H5Dclose( dataset_id );

    H5Gclose( group_id );

    // ========================== Write all species ============================
    for(int iSpec = 0; iSpec < vecSpecies.size(); iSpec++)
    {
        s1 = vecSpecies[iSpec];
        p1 = &( s1->particles );
        string group_name = "/" + s1->species_param.species_type;
        group_id = H5Gcreate(restore_file_id_, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // write the variable "PtclNumber"
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

        // write bmin of species
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

        // write bmax of species
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

        // Write momentum(0,iPart)
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


        // Write momentum(1,iPart)
        dims[0]     = partNumber;
        count[0]    = partNumber;
        offset[0]   = 0;
        stride[0]   = 0;
        block[0]    = 0;

        dataspace_id = H5Screate_simple(RANK, dims, NULL);
        dataset_id = H5Dcreate2(group_id, "momentum1", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(p1->momentum(1,0)) );
        H5Sclose (dataspace_id);
        H5Dclose (dataset_id);


        // Write momentum(2,iPart)
        dims[0]     = partNumber;
        count[0]    = partNumber;
        offset[0]   = 0;
        stride[0]   = 0;
        block[0]    = 0;

        dataspace_id = H5Screate_simple(RANK, dims, NULL);
        dataset_id = H5Dcreate2(group_id, "momentum2", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(p1->momentum(2,0)) );
        H5Sclose (dataspace_id);
        H5Dclose (dataset_id);


        // Write position(0,iPart)
        dims[0]     = partNumber;
        count[0]    = partNumber;
        offset[0]   = 0;
        stride[0]   = 0;
        block[0]    = 0;

        dataspace_id = H5Screate_simple(RANK, dims, NULL);
        dataset_id = H5Dcreate2(group_id, "position0", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(p1->position(0,0)) );
        H5Sclose (dataspace_id);
        H5Dclose (dataset_id);

        if(params.geometry == "2d3v")
        {
            // Write position(1,iPart)
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
void SmileiIO::reloadP( PicParams& params, SmileiMPI* smpi, vector<Species*>& vecSpecies )
{
    Species     *s1;
    Particles   *p1;
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
    int RANK;
    double data_temp;
    hsize_t dims[1];
    hsize_t dims_extented[1];
    hsize_t maxdims[1];
    hsize_t chunk_dims[1];

    RANK = 1;
    maxdims[0] = H5S_UNLIMITED;
    chunk_dims[0] = 2;
    restart = 0;
    data_temp = 66.0;


    // Create the directory "restore" if it not exists
    string dirName = restore_file_dir_;
    if( access(dirName.c_str(), 0) ==-1 )
    {
        mkdir(dirName.c_str(), 0777);
    }


    long long mpi_rk = smpi->getRank();
    string fileName = restore_file_dir_ + "/Restore" + to_string( mpi_rk ) + "_global.h5";
    if( access(fileName.c_str(), 0) ==-1 )  //  The file not exists
    {
        return;
    }

    // Restore000_global.h5
    restore_file_id_  = H5Fopen( fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if(restore_file_id_ < 0)
    {
        WARNING("Can not open restore file: restore_file_id_ = "<<restore_file_id_);
        return;
    }


    // ======================= Read info about the restart ====================
    group_id = H5Gopen(restore_file_id_, "/info", H5P_DEFAULT);
    // Read the variable "timestep"
    dataset_id = H5Dopen(group_id, "timestep", H5P_DEFAULT);
    status = H5Dread (dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &stepStart);
    H5Dclose( dataset_id );

    // Read the variable "restart"
    dataset_id = H5Dopen(group_id, "restart", H5P_DEFAULT);
    status = H5Dread (dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &restart);
    H5Dclose (dataset_id);
    if(restart == 0 && params.is_continue == 0)
    {
        H5Gclose(group_id);
        H5Fclose(restore_file_id_);
        stepStart = 0;
        return;
    }
    else if(restart == 1 && params.is_continue == 0)
    {
        is_restart = 1;
    }
    else if(params.is_continue == 1)
    {
        is_restart = 0;
        stepStart  = 0;
    }


    // Read the variable "nSpecies"
    dataset_id = H5Dopen(group_id, "nSpecies", H5P_DEFAULT);
    status = H5Dread (dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nSpecies);
    H5Dclose( dataset_id );
    H5Gclose( group_id );

    // ====================== Read all stored species ==========================
    for(int iSpec = 0; iSpec < nSpecies; iSpec++)
    {
        s1 = vecSpecies[iSpec];
        p1 = &(s1->particles);
        string group_name = "/" + s1->species_param.species_type;
        group_id = H5Gopen(restore_file_id_, group_name.c_str(), H5P_DEFAULT);

        // Read the variable "PtclNumber" =================================================
        dataset_id = H5Dopen(group_id, "partNumber", H5P_DEFAULT);
        status = H5Dread (dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &partNumber);
        H5Dclose (dataset_id);

        // resize the particles with the partNumber
        p1->initialize(partNumber, params);


        // Read bmin ===================================================================
        dataset_id = H5Dopen(group_id, "bmin", H5P_DEFAULT);
        status = H5Dread (dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(s1->bmin[0]) );
        H5Dclose (dataset_id);

        // Read bmax ===================================================================
        dataset_id = H5Dopen(group_id, "bmax", H5P_DEFAULT);
        status = H5Dread (dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(s1->bmax[0]) );
        H5Dclose (dataset_id);

        // Read momentum(0,iPart) ======================================================
        dataset_id = H5Dopen(group_id, "momentum0", H5P_DEFAULT);
        status = H5Dread (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(p1->momentum(0,0)) );
        H5Dclose (dataset_id);


        // Read momentum(1,iPart) ======================================================
        dataset_id = H5Dopen(group_id, "momentum1", H5P_DEFAULT);
        status = H5Dread (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(p1->momentum(1,0)) );
        H5Dclose (dataset_id);



        // Read momentum(2,iPart)  ======================================================
        dataset_id = H5Dopen(group_id, "momentum2", H5P_DEFAULT);
        status = H5Dread (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(p1->momentum(2,0)) );
        H5Dclose (dataset_id);


        // Read position(0,iPart)  ======================================================
        dataset_id = H5Dopen(group_id, "position0", H5P_DEFAULT);
        status = H5Dread (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(p1->position(0,0)) );
        H5Dclose (dataset_id);

        if(params.geometry == "2d3v")
        {
            // Read position(1,iPart) ======================================================
            dataset_id = H5Dopen(group_id, "position1", H5P_DEFAULT);
            status = H5Dread (dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(p1->position(1,0)) );
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
    string fileName = restore_file_dir_ + "/Restore" + to_string( mpi_rk ) + "_global.h5";

    // Restore000_global.h5
    restore_file_id_  = H5Fopen( fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    string group_name = "/info";
    group_id = H5Gopen(restore_file_id_, group_name.c_str(), H5P_DEFAULT);

    // write the variable "restart" =========================================
    dataset_id = H5Dopen(group_id, "restart", H5P_DEFAULT);
    status = H5Dwrite (dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &restart);

    H5Dclose (dataset_id);
    H5Gclose (group_id);
    H5Fclose (restore_file_id_);
}
