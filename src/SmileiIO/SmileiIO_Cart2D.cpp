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
#include "Species.h"

using namespace std;

SmileiIO_Cart2D::SmileiIO_Cart2D( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies, Diagnostic* diag )
: SmileiIO( params, smpi )
{
    Diagnostic2D* diag2D = static_cast<Diagnostic2D*>(diag);
    reloadP(params, smpi, vecSpecies);
    if(smpi->isMaster())
    {
        if(!is_restart)
        {
            global_file_id_  = H5Fcreate( global_file_name_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        }
        else
        {
            global_file_id_  = H5Fopen( global_file_name_.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        }
        createFieldsPattern(params, smpi, fields);
        createPartsPattern(params, smpi, fields, vecSpecies);
        status = H5Fclose(global_file_id_);
    }

}

SmileiIO_Cart2D::~SmileiIO_Cart2D()
{
}

//> create hdf5 file, datespace, dateset and so on
void SmileiIO_Cart2D::createFieldsPattern( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields )
{

    fieldsGroup.dims_global[3] = params.n_space_global[1] + 1;
    fieldsGroup.dims_global[2] = params.n_space_global[0] + 1;
    fieldsGroup.dims_global[1] = 1;
    fieldsGroup.dims_global[0] = params.n_time / params.dump_step;

    fieldsGroup.ndims_[0] = fieldsGroup.dims_global[0];
    fieldsGroup.ndims_[1] = fieldsGroup.dims_global[1];
    fieldsGroup.ndims_[2] = fieldsGroup.dims_global[2];
    fieldsGroup.ndims_[3] = fieldsGroup.dims_global[3];

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

    // For attribute
    fieldsGroup.aDims = 4;

    createFieldsGroup(fields);


} // END createPattern




void SmileiIO_Cart2D::createPartsPattern( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies )
{
    // For particles, size ofdims_global should be 5: dims_global[nx][ny][nz][nvelocity][ntime]
    // But to be simple, the size is set 4, nz dimension is deleted.
    ptclsGroup.dims_global[3] = vx_dim;
    ptclsGroup.dims_global[2] = params.n_space_global[0];
    ptclsGroup.dims_global[1] = params.n_space_global[1];
    ptclsGroup.dims_global[0] = params.n_time / params.dump_step;

    ptclsGroup.ndims_[0] = ptclsGroup.dims_global[0];
    ptclsGroup.ndims_[1] = ptclsGroup.dims_global[1];
    ptclsGroup.ndims_[2] = ptclsGroup.dims_global[2];
    ptclsGroup.ndims_[3] = ptclsGroup.dims_global[3];

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

    ptclsGroup.aDims = 4;

    createPartsGroup(vecSpecies);

}



void SmileiIO_Cart2D::initVDF( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies )
{
    vx_dim = 200;

    vector<unsigned int> dims_VDF;
    vector<unsigned int> dims_VDF_global;
    dims_VDF.resize(4);
    dims_VDF_global.resize(4);

    dims_VDF[0] = params.n_space[0];
    dims_VDF[1] = params.n_space[1];
    dims_VDF[2] = 1;
    dims_VDF[3] = vx_dim;

    dims_VDF_global[0] = params.n_space_global[0];
    dims_VDF_global[1] = params.n_space_global[1];
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



void SmileiIO_Cart2D::calVDF( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies)
{
    Species *s;
    Particles *p;
    /*

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
    */
}



//! write potential, rho and so on into hdf5 file every some timesteps
void SmileiIO_Cart2D::write( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies, Diagnostic* diag, int itime)
{
    Diagnostic2D* diag2D = static_cast<Diagnostic2D*>(diag);
    if(params.is_calVDF)
    {
        //calVDF( params, smpi, fields, vecSpecies, itime);
    }

    if( itime % params.dump_step == 0 && smpi->isMaster() )
    {
        if(smpi->isMaster())
        {

            ndims_t = itime / params.dump_step - 1;

            // reopen attribute, dataset, dataspace, group, and so on
            global_file_id_  = H5Fopen( global_file_name_.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

            fieldsGroup.group_id = H5Gopen(global_file_id_, "/Fields", H5P_DEFAULT);
            for(int i = 0; i < fieldsGroup.dataset_name.size(); i++)
            {
                fieldsGroup.dataset_id[i] = H5Dopen( fieldsGroup.group_id, fieldsGroup.dataset_name[i], H5P_DEFAULT);
            }

            ptclsGroup.group_id = H5Gopen(global_file_id_, "/VDF", H5P_DEFAULT);
            for(int i = 0; i < ptclsGroup.dataset_name.size(); i++)
            {
                ptclsGroup.dataset_id[i] = H5Dopen( ptclsGroup.group_id, ptclsGroup.dataset_name[i], H5P_DEFAULT);
            }

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



            // write particle velocity distribution function
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

            status = H5Fclose(global_file_id_);
        }
    }



} // END write
