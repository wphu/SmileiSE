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
#include "Species.h"


using namespace std;

SmileiIO_Cart1D::SmileiIO_Cart1D( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies, Diagnostic* diag  )
: SmileiIO( params, smpi )
{
    Diagnostic1D* diag1D = static_cast<Diagnostic1D*>(diag);
    initVDF(params, smpi, fields, vecSpecies);
    if(smpi->isMaster()) {
        createFieldsPattern(params, smpi, fields);
        createPartsPattern(params, smpi, fields, vecSpecies);
        createDiagsPattern(params, smpi, diag1D );
        status = H5Fclose(global_file_id_);
    }
}

SmileiIO_Cart1D::~SmileiIO_Cart1D()
{
}

//> create hdf5 data hierarchical structure: datespace, dateset and so on
void SmileiIO_Cart1D::createFieldsPattern( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields )
{
    hsize_t     dims;

    fieldsGroup.dims_global[3] = params.n_space_global[0] + 1;
    fieldsGroup.dims_global[2] = 1;
    fieldsGroup.dims_global[1] = 1;
    fieldsGroup.dims_global[0] = params.n_time / params.dump_step;

    fieldsGroup.ndims_[0] = fieldsGroup.dims_global[0];
    fieldsGroup.ndims_[1] = fieldsGroup.dims_global[1];
    fieldsGroup.ndims_[2] = fieldsGroup.dims_global[2];
    fieldsGroup.ndims_[3] = fieldsGroup.dims_global[3];


    fieldsGroup.group_id = H5Gcreate(global_file_id_, "/Fields", H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);

    /* Create a datagroup attribute. */
    dims = 4;
    //> dataspace is to descript the structure of data: the number of data dimension and the size of each dimension
    //> the first parameter rank=1: is the number of dimensions used in the dataspace
    fieldsGroup.dataspace_id = H5Screate_simple(1, &dims, NULL);
    fieldsGroup.attribute_id = H5Acreate2 (fieldsGroup.group_id, "dims_global", H5T_STD_I32BE, fieldsGroup.dataspace_id,
                                 H5P_DEFAULT, H5P_DEFAULT);
    /* Write the attribute data. */
    fieldsGroup.status = H5Awrite(fieldsGroup.attribute_id, H5T_NATIVE_INT, fieldsGroup.ndims_);


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
        addField(fields->T_s_global_avg[i]);

    }

    //> if without below process, the method write() will go wrong, no ideas now!!!
    //> output initial 1d_global data===========================================
    data_ =  (double*)malloc(fieldsGroup.dims_global[3] * fieldsGroup.dims_global[2] * fieldsGroup.dims_global[1] * fieldsGroup.dims_global[0] * sizeof(double));
    for( int i = 0; i < fieldsGroup.dims_global[3] * fieldsGroup.dims_global[2] * fieldsGroup.dims_global[1] * fieldsGroup.dims_global[0]; i++)
    {
      data_[i] = 20.0;
    }

    for(int i = 0; i < fieldsGroup.dataset_id.size(); i++)
    {
        fieldsGroup.status = H5Dwrite(fieldsGroup.dataset_id[i], H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_);
    }
    free(data_);


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



    // close attribute, dataset, dataspace, group, and so on
    fieldsGroup.status = H5Aclose( fieldsGroup.attribute_id );
    fieldsGroup.status = H5Sclose( fieldsGroup.dataspace_id );
    for(int i = 0; i < fieldsGroup.dataset_id.size(); i++)
    {
        fieldsGroup.status = H5Dclose( fieldsGroup.dataset_id[i] );
    }
    fieldsGroup.status = H5Gclose( fieldsGroup.group_id );


} // END createPattern


void SmileiIO_Cart1D::createPartsPattern( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies )
{
    hsize_t     dims;

    // For particles, size ofdims_global should be 5: dims_global[nx][ny][nz][nvelocity][ntime]
    // But to be simple, the size is set 4, nz dimension is deleted.
    ptclsGroup.dims_global[3] = vx_dim;
    ptclsGroup.dims_global[2] = params.n_space_global[0];
    ptclsGroup.dims_global[1] = 1;
    ptclsGroup.dims_global[0] = params.n_time / params.dump_step;

    ptclsGroup.ndims_[0] = ptclsGroup.dims_global[0];
    ptclsGroup.ndims_[1] = ptclsGroup.dims_global[1];
    ptclsGroup.ndims_[2] = ptclsGroup.dims_global[2];
    ptclsGroup.ndims_[3] = ptclsGroup.dims_global[3];


    ptclsGroup.group_id = H5Gcreate(global_file_id_, "/VDF", H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);

    /* Create a datagroup attribute. */
    dims = 4;
    //> dataspace is to descript the structure of data: the number of data dimension and the size of each dimension
    //> the first parameter rank=1: is the number of dimensions used in the dataspace
    ptclsGroup.dataspace_id = H5Screate_simple(1, &dims, NULL);
    ptclsGroup.attribute_id = H5Acreate2 (ptclsGroup.group_id, "dims_global", H5T_STD_I32BE, ptclsGroup.dataspace_id,
                                 H5P_DEFAULT, H5P_DEFAULT);
    /* Write the attribute data. */
    ptclsGroup.status = H5Awrite(ptclsGroup.attribute_id, H5T_NATIVE_INT, ptclsGroup.ndims_);



    Species *s;
    Particles *p;


    for(int isp=0; isp<vx_VDF.size(); isp++)
    {
        s = vecSpecies[isp];
        p = &(s->particles);

        string fieldName = "VDF_" + s->species_param.species_type;
        ptclsGroup.dataset_stringName.push_back(fieldName);
        const char* name = ptclsGroup.dataset_stringName.back().c_str();
        ptclsGroup.dataset_name.push_back(name);

        /* Create the data space for the dataset. */
        ptclsGroup.dataspace_id = H5Screate_simple(4, ptclsGroup.dims_global, NULL);
        hid_t id = H5Dcreate2(ptclsGroup.group_id, ptclsGroup.dataset_name.back(), H5T_NATIVE_DOUBLE, ptclsGroup.dataspace_id,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        ptclsGroup.dataset_id.push_back(id);

        ptclsGroup.dataset_data.push_back(vx_VDF_global[isp]->data_);

        // for VDF_tot_ the dimension is same with VDF_, but only [nx=0] is meaningful
        fieldName = "VDF_tot_" + s->species_param.species_type;
        ptclsGroup.dataset_stringName.push_back(fieldName);
        name = ptclsGroup.dataset_stringName.back().c_str();
        ptclsGroup.dataset_name.push_back(name);

        /* Create the data space for the dataset. */
        ptclsGroup.dataspace_id = H5Screate_simple(4, ptclsGroup.dims_global, NULL);
        id = H5Dcreate2(ptclsGroup.group_id, ptclsGroup.dataset_name.back(), H5T_NATIVE_DOUBLE, ptclsGroup.dataspace_id,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        ptclsGroup.dataset_id.push_back(id);

        ptclsGroup.dataset_data.push_back(vx_VDF_tot_global[isp]->data_);


    }


    //> if without below process, the method write() will go wrong, no ideas now!!!
    //> output initial 1d_global data===========================================
    data_ =  (double*)malloc(ptclsGroup.dims_global[3] * ptclsGroup.dims_global[2] * ptclsGroup.dims_global[1] * ptclsGroup.dims_global[0] * sizeof(double));
    for( int i = 0; i < ptclsGroup.dims_global[3] * ptclsGroup.dims_global[2] * ptclsGroup.dims_global[1] * ptclsGroup.dims_global[0]; i++)
    {
      data_[i] = 20.0;
    }

    for(int i = 0; i < ptclsGroup.dataset_id.size(); i++)
    {
        ptclsGroup.status = H5Dwrite(ptclsGroup.dataset_id[i], H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_);
    }
    free(data_);


    ptclsGroup.offset[0] = 0;
    ptclsGroup.offset[1] = 0;
    ptclsGroup.offset[2] = 0;
    ptclsGroup.offset[3] = 0;

    ptclsGroup.stride[0] = 1;
    ptclsGroup.stride[1] = 1;
    ptclsGroup.stride[2] = 1;
    ptclsGroup.stride[3] = 1;

    ptclsGroup.block[0] = 1;
    ptclsGroup.block[1] = 1;
    ptclsGroup.block[2] = 1;
    ptclsGroup.block[3] = 1;


    // close attribute, dataset, dataspace, group, and so on
    ptclsGroup.status = H5Aclose( ptclsGroup.attribute_id );
    ptclsGroup.status = H5Sclose( ptclsGroup.dataspace_id );
    for(int i = 0; i < ptclsGroup.dataset_id.size(); i++)
    {
        ptclsGroup.status = H5Dclose( ptclsGroup.dataset_id[i] );
    }
    ptclsGroup.status = H5Gclose( ptclsGroup.group_id );

}


// Create particles h5 file pattern
void SmileiIO_Cart1D::createDiagsPattern( PicParams& params, SmileiMPI* smpi, Diagnostic1D* diag1D )
{
    string diagName;
    const char* name;
    hid_t dataset_id;

    diagsGroup.group_id = H5Gcreate(global_file_id_, "/Diagnostic", H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);

    // =======create dataset for particleFlux================================
    diagsGroup.dims_global[3] = 2;
    diagsGroup.dims_global[2] = diag1D->n_species;
    diagsGroup.dims_global[1] = 1;
    diagsGroup.dims_global[0] = params.n_time / params.dump_step;

    // dataset name
    diagName = "particleFlux";
    diagsGroup.dataset_stringName.push_back(diagName);
    name = diagsGroup.dataset_stringName.back().c_str();
    diagsGroup.dataset_name.push_back(name);

    // Create the data space for the dataset
    diagsGroup.dataspace_id = H5Screate_simple(4, diagsGroup.dims_global, NULL);
    // Create the dataset
    dataset_id = H5Dcreate2(diagsGroup.group_id, diagsGroup.dataset_name.back(), H5T_NATIVE_DOUBLE, diagsGroup.dataspace_id,
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    diagsGroup.dataset_id.push_back(dataset_id);

    // Write initial data to the dataset
    data_ =  (double*)malloc(diagsGroup.dims_global[3] * diagsGroup.dims_global[2] * diagsGroup.dims_global[1] * diagsGroup.dims_global[0] * sizeof(double));
    for( int i = 0; i < diagsGroup.dims_global[3] * diagsGroup.dims_global[2] * diagsGroup.dims_global[1] * diagsGroup.dims_global[0]; i++)
    {
      data_[i] = 2.0;
    }
    diagsGroup.status = H5Dwrite(diagsGroup.dataset_id.back(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_);
    free(data_);


    // =======create dataset for heatFlux================================
    diagsGroup.dims_global[3] = 2;
    diagsGroup.dims_global[2] = diag1D->n_species;
    diagsGroup.dims_global[1] = 1;
    diagsGroup.dims_global[0] = params.n_time / params.dump_step;

    // dataset name
    diagName = "heatFlux";
    diagsGroup.dataset_stringName.push_back(diagName);
    name = diagsGroup.dataset_stringName.back().c_str();
    diagsGroup.dataset_name.push_back(name);

    // Create the data space for the dataset
    diagsGroup.dataspace_id = H5Screate_simple(4, diagsGroup.dims_global, NULL);
    // Create the dataset
    dataset_id = H5Dcreate2(diagsGroup.group_id, diagsGroup.dataset_name.back(), H5T_NATIVE_DOUBLE, diagsGroup.dataspace_id,
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    diagsGroup.dataset_id.push_back(dataset_id);

    // Write initial data to the dataset
    data_ =  (double*)malloc(diagsGroup.dims_global[3] * diagsGroup.dims_global[2] * diagsGroup.dims_global[1] * diagsGroup.dims_global[0] * sizeof(double));
    for( int i = 0; i < diagsGroup.dims_global[3] * diagsGroup.dims_global[2] * diagsGroup.dims_global[1] * diagsGroup.dims_global[0]; i++)
    {
      data_[i] = 2.0;
    }
    diagsGroup.status = H5Dwrite(diagsGroup.dataset_id.back(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_);
    free(data_);


    // =======create dataset for angleDist================================
    diagsGroup.dims_global[3] = 90;
    diagsGroup.dims_global[2] = 2;
    diagsGroup.dims_global[1] = diag1D->n_species;
    diagsGroup.dims_global[0] = params.n_time / params.dump_step;

    // dataset name
    diagName = "angleDist";
    diagsGroup.dataset_stringName.push_back(diagName);
    name = diagsGroup.dataset_stringName.back().c_str();
    diagsGroup.dataset_name.push_back(name);

    // Create the data space for the dataset
    diagsGroup.dataspace_id = H5Screate_simple(4, diagsGroup.dims_global, NULL);
    // Create the dataset
    dataset_id = H5Dcreate2(diagsGroup.group_id, diagsGroup.dataset_name.back(), H5T_NATIVE_DOUBLE, diagsGroup.dataspace_id,
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    diagsGroup.dataset_id.push_back(dataset_id);

    // Write initial data to the dataset
    data_ =  (double*)malloc(diagsGroup.dims_global[3] * diagsGroup.dims_global[2] * diagsGroup.dims_global[1] * diagsGroup.dims_global[0] * sizeof(double));
    for( int i = 0; i < diagsGroup.dims_global[3] * diagsGroup.dims_global[2] * diagsGroup.dims_global[1] * diagsGroup.dims_global[0]; i++)
    {
      data_[i] = 2.0;
    }
    diagsGroup.status = H5Dwrite(diagsGroup.dataset_id.back(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_);
    free(data_);

    diagsGroup.stride[0] = 1;
    diagsGroup.stride[1] = 1;
    diagsGroup.stride[2] = 1;
    diagsGroup.stride[3] = 1;

    diagsGroup.block[0] = 1;
    diagsGroup.block[1] = 1;
    diagsGroup.block[2] = 1;
    diagsGroup.block[3] = 1;

    // close dataset and group
    for(int i = 0; i < diagsGroup.dataset_id.size(); i++)
    {
        diagsGroup.status = H5Dclose( diagsGroup.dataset_id[i] );
    }
    diagsGroup.status = H5Gclose( diagsGroup.group_id );
}



void SmileiIO_Cart1D::initVDF( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies )
{
    vx_dim = 200;

    vector<unsigned int> dims_VDF;
    vector<unsigned int> dims_VDF_global;
    dims_VDF.resize(4);
    dims_VDF_global.resize(4);

    dims_VDF[0] = params.n_space[0];
    dims_VDF[1] = 1;
    dims_VDF[2] = 1;
    dims_VDF[3] = vx_dim;

    dims_VDF_global[0] = params.n_space_global[0];
    dims_VDF_global[1] = 1;
    dims_VDF_global[2] = 1;
    dims_VDF_global[3] = vx_dim;

    for(int isp=0; isp<vecSpecies.size(); isp++)
    {
        vx_VDF.push_back(new Array4D());
        vx_VDF[isp]->allocateDims(dims_VDF);

        vx_VDF_global.push_back(new Array4D());
        vx_VDF_global[isp]->allocateDims(dims_VDF_global);

        vx_VDF_tot_global.push_back(new Array4D());
        vx_VDF_tot_global[isp]->allocateDims(dims_VDF_global);
    }

}



void SmileiIO_Cart1D::calVDF( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies)
{
    Species *s;
    Particles *p;


    for(int isp=0; isp<vx_VDF.size(); isp++)
    {
        s = vecSpecies[isp];
        p = &(s->particles);

        vector<double> x_cell(3,0);
        x_cell[0] = 0;
        x_cell[1] = 0;
        x_cell[2] = 0;
        vxMax = 3*sqrt(2.0 * s->species_param.thermT[0] * params.const_e / s->species_param.mass);
        //WARNING("thermalVelocity" <<  s->species_param.thermT[0] );
        vxMin = -vxMax;
        vx_d = (vxMax - vxMin) / vx_dim;
        int vx_dim2 = vx_dim / 2;

        vx_VDF[isp]->put_to(0.0);
        for(int ibin = 0; ibin < ( s->bmin.size() ); ibin++)
        {
            for(int iPart = s->bmin[ibin]; iPart < s->bmax[ibin]; iPart++)
            {
                int ivx = p->momentum(0,iPart) / vx_d + vx_dim2;
                if( ivx < 0 ) {
                    ivx = 0;
                }
                if( ivx >= vx_dim ) {
                    ivx = vx_dim - 1;
                }
                (*vx_VDF[isp])(ibin,0,0,ivx) += 1.0;
            }
        }
        smpi->gatherVDF(vx_VDF_global[isp], vx_VDF[isp]);

        vx_VDF_tot_global[isp]->put_to(0.0);
        for (int ibin = 0; ibin < vx_VDF_global[isp]->dims_[0]; ibin++)
        {
            for (int ivx = 0; ivx < vx_VDF_global[isp]->dims_[3]; ivx++)
            {
                (*vx_VDF_tot_global[isp])(0,0,0,ivx) += (*vx_VDF_global[isp])(ibin,0,0,ivx);
            }

        }


    }
}



//! write potential, rho and so on into hdf5 file every some timesteps
void SmileiIO_Cart1D::write( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies, Diagnostic* diag)
{

    Diagnostic1D* diag1D = static_cast<Diagnostic1D*>(diag);
    calVDF( params, smpi, fields, vecSpecies);

    if(smpi->isMaster()) {

        // reopen attribute, dataset, dataspace, group, and so on
        global_file_id_  = H5Fopen( "Fields_global.h5", H5F_ACC_RDWR, H5P_DEFAULT);


        // =============write fields============================================
        fieldsGroup.group_id = H5Gopen(global_file_id_, "/Fields", H5P_DEFAULT);
        for(int i = 0; i < fieldsGroup.dataset_name.size(); i++)
        {
            fieldsGroup.dataset_id[i] = H5Dopen( fieldsGroup.group_id, fieldsGroup.dataset_name[i], H5P_DEFAULT);
        }

        fieldsGroup.offset[0] = ndims_t;
        fieldsGroup.count[0]  = 1;
        fieldsGroup.count[1]  = fieldsGroup.dims_global[1];
        fieldsGroup.count[2]  = fieldsGroup.dims_global[2];
        fieldsGroup.count[3]  = fieldsGroup.dims_global[3];

        for(int i = 0; i < fieldsGroup.dataset_id.size(); i++)
        {
            fieldsGroup.memspace_id = H5Screate_simple (4, fieldsGroup.count, NULL);
            fieldsGroup.dataspace_id = H5Dget_space (fieldsGroup.dataset_id[i]);
            // H5Sselect_hyperslab: define the selection of subset in the dataspace
            fieldsGroup.status = H5Sselect_hyperslab (fieldsGroup.dataspace_id, H5S_SELECT_SET, fieldsGroup.offset,
                                              fieldsGroup.stride, fieldsGroup.count, fieldsGroup.block);
            fieldsGroup.status = H5Dwrite (fieldsGroup.dataset_id[i], H5T_NATIVE_DOUBLE, fieldsGroup.memspace_id,
                                 fieldsGroup.dataspace_id, H5P_DEFAULT, fieldsGroup.dataset_data[i]);

            fieldsGroup.status = H5Sclose (fieldsGroup.memspace_id);
            fieldsGroup.status = H5Sclose (fieldsGroup.dataspace_id);
        }

        // close dataset, group, and so on
        for(int i = 0; i < fieldsGroup.dataset_id.size(); i++)
        {
            fieldsGroup.status = H5Dclose( fieldsGroup.dataset_id[i] );
        }
        fieldsGroup.status = H5Gclose( fieldsGroup.group_id );



        // ==============write particle velocity distribution function=========================
        ptclsGroup.group_id = H5Gopen(global_file_id_, "/VDF", H5P_DEFAULT);
        for(int i = 0; i < ptclsGroup.dataset_name.size(); i++)
        {
            ptclsGroup.dataset_id[i] = H5Dopen( ptclsGroup.group_id, ptclsGroup.dataset_name[i], H5P_DEFAULT);
        }

        ptclsGroup.offset[0] = ndims_t;
        ptclsGroup.count[0]  = 1;
        ptclsGroup.count[1]  = ptclsGroup.dims_global[1];
        ptclsGroup.count[2]  = ptclsGroup.dims_global[2];
        ptclsGroup.count[3]  = ptclsGroup.dims_global[3];
        for(int i = 0; i < ptclsGroup.dataset_id.size(); i++)
        {
            ptclsGroup.memspace_id = H5Screate_simple (4, ptclsGroup.count, NULL);
            ptclsGroup.dataspace_id = H5Dget_space (ptclsGroup.dataset_id[i]);
            ptclsGroup.status = H5Sselect_hyperslab (ptclsGroup.dataspace_id, H5S_SELECT_SET, ptclsGroup.offset,
                                              ptclsGroup.stride, ptclsGroup.count, ptclsGroup.block);
            ptclsGroup.status = H5Dwrite (ptclsGroup.dataset_id[i], H5T_NATIVE_DOUBLE, ptclsGroup.memspace_id,
                                 ptclsGroup.dataspace_id, H5P_DEFAULT, ptclsGroup.dataset_data[i]);

            ptclsGroup.status = H5Sclose (ptclsGroup.memspace_id);
            ptclsGroup.status = H5Sclose (ptclsGroup.dataspace_id);
        }

        // close dataset, group, and so on
        for(int i = 0; i < ptclsGroup.dataset_id.size(); i++)
        {
            ptclsGroup.status = H5Dclose( ptclsGroup.dataset_id[i] );
        }
        ptclsGroup.status = H5Gclose( ptclsGroup.group_id );


        // ==============write Diagnostic: particleFlux, heatFlux and angleDist============
        diagsGroup.group_id = H5Gopen(global_file_id_, "/Diagnostic", H5P_DEFAULT);
        for(int i = 0; i < diagsGroup.dataset_name.size(); i++)
        {
            diagsGroup.dataset_id[i] = H5Dopen( diagsGroup.group_id, diagsGroup.dataset_name[i], H5P_DEFAULT);
        }

        for(int ispec = 0; ispec < diag1D->n_species; ispec++)
        {
            for(int iDirection = 0; iDirection < 2; iDirection++)
            {
                diagsGroup.offset[0] = ndims_t;
                diagsGroup.offset[1] = 0;
                diagsGroup.offset[2] = ispec;
                diagsGroup.offset[3] = iDirection;
                diagsGroup.count[0]  = 1;
                diagsGroup.count[1]  = 1;
                diagsGroup.count[2]  = 1;
                diagsGroup.count[3]  = 1;

                // particleFlux
                diagsGroup.memspace_id = H5Screate_simple (4, diagsGroup.count, NULL);
                diagsGroup.dataspace_id = H5Dget_space (diagsGroup.dataset_id[0]);
                diagsGroup.status = H5Sselect_hyperslab (diagsGroup.dataspace_id, H5S_SELECT_SET, diagsGroup.offset,
                                                  diagsGroup.stride, diagsGroup.count, diagsGroup.block);
                diagsGroup.status = H5Dwrite (diagsGroup.dataset_id[0], H5T_NATIVE_DOUBLE, diagsGroup.memspace_id,
                                     diagsGroup.dataspace_id, H5P_DEFAULT, &(diag1D->particleFlux[ispec][iDirection]) );

                diagsGroup.status = H5Sclose (diagsGroup.memspace_id);
                diagsGroup.status = H5Sclose (diagsGroup.dataspace_id);

                // heatFlux
                diagsGroup.memspace_id = H5Screate_simple (4, diagsGroup.count, NULL);
                diagsGroup.dataspace_id = H5Dget_space (diagsGroup.dataset_id[1]);
                diagsGroup.status = H5Sselect_hyperslab (diagsGroup.dataspace_id, H5S_SELECT_SET, diagsGroup.offset,
                                                  diagsGroup.stride, diagsGroup.count, diagsGroup.block);
                diagsGroup.status = H5Dwrite (diagsGroup.dataset_id[1], H5T_NATIVE_DOUBLE, diagsGroup.memspace_id,
                                     diagsGroup.dataspace_id, H5P_DEFAULT, &(diag1D->heatFlux[ispec][iDirection]) );

                diagsGroup.status = H5Sclose (diagsGroup.memspace_id);
                diagsGroup.status = H5Sclose (diagsGroup.dataspace_id);


                diagsGroup.offset[0] = ndims_t;
                diagsGroup.offset[1] = ispec;
                diagsGroup.offset[2] = iDirection;
                diagsGroup.offset[3] = 0;
                diagsGroup.count[0]  = 1;
                diagsGroup.count[1]  = 1;
                diagsGroup.count[2]  = 1;
                diagsGroup.count[3]  = 90;

                // angleDist
                diagsGroup.memspace_id = H5Screate_simple (4, diagsGroup.count, NULL);
                diagsGroup.dataspace_id = H5Dget_space (diagsGroup.dataset_id[2]);
                diagsGroup.status = H5Sselect_hyperslab (diagsGroup.dataspace_id, H5S_SELECT_SET, diagsGroup.offset,
                                                  diagsGroup.stride, diagsGroup.count, diagsGroup.block);
                diagsGroup.status = H5Dwrite (diagsGroup.dataset_id[2], H5T_NATIVE_DOUBLE, diagsGroup.memspace_id,
                                     diagsGroup.dataspace_id, H5P_DEFAULT, &(diag1D->angleDist[ispec][iDirection][0]) );

                diagsGroup.status = H5Sclose (diagsGroup.memspace_id);
                diagsGroup.status = H5Sclose (diagsGroup.dataspace_id);

            }
        }


        // close dataset, group, and so on
        for(int i = 0; i < diagsGroup.dataset_id.size(); i++)
        {
            diagsGroup.status = H5Dclose( diagsGroup.dataset_id[i] );
        }
        diagsGroup.status = H5Gclose( diagsGroup.group_id );


        ndims_t++;
        status = H5Fclose(global_file_id_);
    }


} // END write
