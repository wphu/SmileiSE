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
    reloadP(params, smpi, vecSpecies);
    if(smpi->isMaster()) {
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
    fieldsGroup.dims_global[3] = params.n_space_global[0] + 1;
    fieldsGroup.dims_global[2] = 1;
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


void SmileiIO_Cart1D::createPartsPattern( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies )
{

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


// Create particles h5 file pattern
void SmileiIO_Cart1D::createDiagsPattern( PicParams& params, SmileiMPI* smpi, Diagnostic1D* diag1D )
{
    string diagName;
    const char* name;
    hid_t dataset_id;

    // =======set stride and block, and close dataset and group=================
    diagsGroup.stride[0] = 1;
    diagsGroup.stride[1] = 1;
    diagsGroup.stride[2] = 1;
    diagsGroup.stride[3] = 1;

    diagsGroup.block[0] = 1;
    diagsGroup.block[1] = 1;
    diagsGroup.block[2] = 1;
    diagsGroup.block[3] = 1;

    // ======= create diagsGroup ================================
    // dataset name
    diagName = "particleFlux";
    diagsGroup.dataset_stringName.push_back(diagName);
    name = diagsGroup.dataset_stringName.back().c_str();
    diagsGroup.dataset_name.push_back(name);

    // dataset name
    diagName = "heatFlux";
    diagsGroup.dataset_stringName.push_back(diagName);
    name = diagsGroup.dataset_stringName.back().c_str();
    diagsGroup.dataset_name.push_back(name);

    // dataset name
    diagName = "angleDist";
    diagsGroup.dataset_stringName.push_back(diagName);
    name = diagsGroup.dataset_stringName.back().c_str();
    diagsGroup.dataset_name.push_back(name);

    // dataset name
    diagName = "particleNumber";
    diagsGroup.dataset_stringName.push_back(diagName);
    name = diagsGroup.dataset_stringName.back().c_str();
    diagsGroup.dataset_name.push_back(name);

    // dataset name
    diagName = "kineticEnergy";
    diagsGroup.dataset_stringName.push_back(diagName);
    name = diagsGroup.dataset_stringName.back().c_str();
    diagsGroup.dataset_name.push_back(name);


    if(!is_restart)
    {
        diagsGroup.group_id = H5Gcreate(global_file_id_, "/Diagnostic", H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);
        // ======= Create hdf5 struct for "/Diagnostic" group ================================
        int iDiag = 0;
        // Create dataset for particleFlux
        diagsGroup.dims_global[3] = 2;
        diagsGroup.dims_global[2] = diag1D->n_species;
        diagsGroup.dims_global[1] = 1;
        diagsGroup.dims_global[0] = params.n_time / params.dump_step;

        diagsGroup.dataspace_id = H5Screate_simple(4, diagsGroup.dims_global, NULL);
        dataset_id = H5Dcreate2(diagsGroup.group_id, diagsGroup.dataset_name[iDiag], H5T_NATIVE_DOUBLE, diagsGroup.dataspace_id,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        diagsGroup.dataset_id.push_back(dataset_id);
        H5Dclose( diagsGroup.dataset_id[iDiag] );
        H5Sclose( diagsGroup.dataspace_id );


        // Create dataset for heatFlux
        diagsGroup.dims_global[3] = 2;
        diagsGroup.dims_global[2] = diag1D->n_species;
        diagsGroup.dims_global[1] = 1;
        diagsGroup.dims_global[0] = params.n_time / params.dump_step;

        iDiag++;
        diagsGroup.dataspace_id = H5Screate_simple(4, diagsGroup.dims_global, NULL);
        dataset_id = H5Dcreate2(diagsGroup.group_id, diagsGroup.dataset_name[iDiag], H5T_NATIVE_DOUBLE, diagsGroup.dataspace_id,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        diagsGroup.dataset_id.push_back(dataset_id);
        H5Dclose( diagsGroup.dataset_id[iDiag] );
        H5Sclose( diagsGroup.dataspace_id );



        // Create dataset for angleDist
        diagsGroup.dims_global[3] = 90;
        diagsGroup.dims_global[2] = 2;
        diagsGroup.dims_global[1] = diag1D->n_species;
        diagsGroup.dims_global[0] = params.n_time / params.dump_step;

        iDiag++;
        diagsGroup.dataspace_id = H5Screate_simple(4, diagsGroup.dims_global, NULL);
        dataset_id = H5Dcreate2(diagsGroup.group_id, diagsGroup.dataset_name[iDiag], H5T_NATIVE_DOUBLE, diagsGroup.dataspace_id,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        diagsGroup.dataset_id.push_back(dataset_id);
        H5Dclose( diagsGroup.dataset_id[iDiag] );
        H5Sclose( diagsGroup.dataspace_id );


        // Create dataset for particleNumber
        diagsGroup.dims_global[3] = diag1D->n_species;
        diagsGroup.dims_global[2] = 1;
        diagsGroup.dims_global[1] = 1;
        diagsGroup.dims_global[0] = params.n_time / params.dump_step;

        iDiag++;
        diagsGroup.dataspace_id = H5Screate_simple(4, diagsGroup.dims_global, NULL);
        dataset_id = H5Dcreate2(diagsGroup.group_id, diagsGroup.dataset_name[iDiag], H5T_NATIVE_INT, diagsGroup.dataspace_id,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        diagsGroup.dataset_id.push_back(dataset_id);
        H5Dclose( diagsGroup.dataset_id[iDiag] );
        H5Sclose( diagsGroup.dataspace_id );


        // Create dataset for particleNumber
        diagsGroup.dims_global[3] = diag1D->n_species;
        diagsGroup.dims_global[2] = 1;
        diagsGroup.dims_global[1] = 1;
        diagsGroup.dims_global[0] = params.n_time / params.dump_step;

        iDiag++;
        diagsGroup.dataspace_id = H5Screate_simple(4, diagsGroup.dims_global, NULL);
        dataset_id = H5Dcreate2(diagsGroup.group_id, diagsGroup.dataset_name[iDiag], H5T_NATIVE_DOUBLE, diagsGroup.dataspace_id,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        diagsGroup.dataset_id.push_back(dataset_id);
        H5Dclose( diagsGroup.dataset_id[iDiag] );
        H5Sclose( diagsGroup.dataspace_id );

        diagsGroup.status = H5Gclose( diagsGroup.group_id );
    }
    diagsGroup.dataset_id.resize( diagsGroup.dataset_name.size() );
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



void SmileiIO_Cart1D::calVDF( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies, int itime)
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

        if( (itime % (params.dump_step + 1)) == 0 )
        {
            vx_VDF_tot_global[isp]->put_to(0.0);
        }

        for (int ibin = 0; ibin < vx_VDF_global[isp]->dims_[0]; ibin++)
        {
            for (int ivx = 0; ivx < vx_VDF_global[isp]->dims_[3]; ivx++)
            {
                (*vx_VDF_tot_global[isp])(0,0,0,ivx) += (*vx_VDF_global[isp])(ibin,0,0,ivx);
            }

        }

        if( (itime % (params.dump_step)) == 0 )
        {
            for (int ivx = 0; ivx < vx_VDF_global[isp]->dims_[3]; ivx++)
            {
                (*vx_VDF_tot_global[isp])(0,0,0,ivx) /= params.dump_step;
            }
        }

    }
}



//! write potential, rho and so on into hdf5 file every some timesteps
void SmileiIO_Cart1D::write( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies, Diagnostic* diag, int itime)
{

    Diagnostic1D* diag1D = static_cast<Diagnostic1D*>(diag);
    calVDF( params, smpi, fields, vecSpecies, itime);

    if(smpi->isMaster()) {

        ndims_t = itime / params.dump_step - 1;

        // reopen attribute, dataset, dataspace, group, and so on
        global_file_id_  = H5Fopen( global_file_name_.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

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
        // ==============particleNumber, kineticEnergy                         ============
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

        diagsGroup.offset[0] = ndims_t;
        diagsGroup.offset[1] = 0;
        diagsGroup.offset[2] = 0;
        diagsGroup.offset[3] = 0;
        diagsGroup.count[0]  = 1;
        diagsGroup.count[1]  = 1;
        diagsGroup.count[2]  = 1;
        diagsGroup.count[3]  = diag1D->n_species;

        // particleNumber
        diagsGroup.memspace_id = H5Screate_simple (4, diagsGroup.count, NULL);
        diagsGroup.dataspace_id = H5Dget_space (diagsGroup.dataset_id[3]);
        diagsGroup.status = H5Sselect_hyperslab (diagsGroup.dataspace_id, H5S_SELECT_SET, diagsGroup.offset,
                                          diagsGroup.stride, diagsGroup.count, diagsGroup.block);
        diagsGroup.status = H5Dwrite (diagsGroup.dataset_id[3], H5T_NATIVE_INT, diagsGroup.memspace_id,
                             diagsGroup.dataspace_id, H5P_DEFAULT, &(diag1D->particleNumber[0]) );

        diagsGroup.status = H5Sclose (diagsGroup.memspace_id);
        diagsGroup.status = H5Sclose (diagsGroup.dataspace_id);


        // kineticEnergy
        diagsGroup.memspace_id = H5Screate_simple (4, diagsGroup.count, NULL);
        diagsGroup.dataspace_id = H5Dget_space (diagsGroup.dataset_id[4]);
        diagsGroup.status = H5Sselect_hyperslab (diagsGroup.dataspace_id, H5S_SELECT_SET, diagsGroup.offset,
                                          diagsGroup.stride, diagsGroup.count, diagsGroup.block);
        diagsGroup.status = H5Dwrite (diagsGroup.dataset_id[4], H5T_NATIVE_DOUBLE, diagsGroup.memspace_id,
                             diagsGroup.dataspace_id, H5P_DEFAULT, &(diag1D->kineticEnergy[0]) );

        diagsGroup.status = H5Sclose (diagsGroup.memspace_id);
        diagsGroup.status = H5Sclose (diagsGroup.dataspace_id);


        // close dataset, group, and so on
        for(int i = 0; i < diagsGroup.dataset_id.size(); i++)
        {
            diagsGroup.status = H5Dclose( diagsGroup.dataset_id[i] );
        }
        diagsGroup.status = H5Gclose( diagsGroup.group_id );

        status = H5Fclose(global_file_id_);
    }


} // END write
