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

    if( smpi->isMaster() )
    {
        // create data patterns
        createFieldsPattern(params, fields);
        // at present, not calculate VDF in the core, VDF can be calculated from particles in the "restore" directory
        //createPartsPattern(params, fields, vecSpecies);

    }

}

SmileiIO_Cart1D::~SmileiIO_Cart1D()
{
}

//> create hdf5 data hierarchical structure: datespace, dateset and so on
void SmileiIO_Cart1D::createFieldsPattern( PicParams& params, ElectroMagn* fields )
{
    fieldsGroup.dims_global[2] = params.n_space_global[0] + 1;
    fieldsGroup.dims_global[1] = 1;
    fieldsGroup.dims_global[0] = 1;

    fieldsGroup.ndims_[0] = fieldsGroup.dims_global[0];
    fieldsGroup.ndims_[1] = fieldsGroup.dims_global[1];
    fieldsGroup.ndims_[2] = fieldsGroup.dims_global[2];


    fieldsGroup.offset[0] = 0;
    fieldsGroup.offset[1] = 0;
    fieldsGroup.offset[2] = 0;


    fieldsGroup.stride[0] = 1;
    fieldsGroup.stride[1] = 1;
    fieldsGroup.stride[2] = 1;


    fieldsGroup.block[0] = 1;
    fieldsGroup.block[1] = 1;
    fieldsGroup.block[2] = 1;


    // For attribute
    fieldsGroup.aDims = 3;

    createFieldsGroup(fields);

} // END createPattern

// Create particles h5 file pattern
void SmileiIO_Cart1D::createDiagsPattern( PicParams& params, Diagnostic1D* diag1D )
{
    string diagName;
    const char* h5_name;
    hid_t dataset_id;
    int data_dims = 3;

    // =======set stride and block, and close dataset and group=================
    diagsGroup.stride[0] = 1;
    diagsGroup.stride[1] = 1;
    diagsGroup.stride[2] = 1;


    diagsGroup.block[0] = 1;
    diagsGroup.block[1] = 1;
    diagsGroup.block[2] = 1;


    // ======= create diagsGroup ================================
    // dataset name
    diagName = "particleFlux";
    diagsGroup.dataset_stringName.push_back(diagName);


    // dataset name
    diagName = "heatFlux";
    diagsGroup.dataset_stringName.push_back(diagName);


    // dataset name
    diagName = "angleDist";
    diagsGroup.dataset_stringName.push_back(diagName);


    // dataset name
    diagName = "particleNumber";
    diagsGroup.dataset_stringName.push_back(diagName);


    // dataset name
    diagName = "kineticEnergy";
    diagsGroup.dataset_stringName.push_back(diagName);



    if(!is_restart)
    {
        diagsGroup.group_id = H5Gcreate(data_file_id, "/Diagnostic", H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);
        // ======= Create hdf5 struct for "/Diagnostic" group ================================
        int iDiag = 0;
        // Create dataset for particleFlux
        diagsGroup.dims_global[2] = 2;
        diagsGroup.dims_global[1] = diag1D->n_species;
        diagsGroup.dims_global[0] = 1;


        diagsGroup.dataspace_id = H5Screate_simple(data_dims, diagsGroup.dims_global, NULL);
        h5_name = diagsGroup.dataset_stringName[iDiag].c_str();
        dataset_id = H5Dcreate2(diagsGroup.group_id, h5_name, H5T_NATIVE_DOUBLE, diagsGroup.dataspace_id,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        diagsGroup.dataset_id.push_back(dataset_id);
        H5Dclose( diagsGroup.dataset_id[iDiag] );
        H5Sclose( diagsGroup.dataspace_id );


        // Create dataset for heatFlux
        diagsGroup.dims_global[2] = 2;
        diagsGroup.dims_global[1] = diag1D->n_species;
        diagsGroup.dims_global[0] = 1;

        iDiag++;
        diagsGroup.dataspace_id = H5Screate_simple(data_dims, diagsGroup.dims_global, NULL);
        h5_name = diagsGroup.dataset_stringName[iDiag].c_str();
        dataset_id = H5Dcreate2(diagsGroup.group_id, h5_name, H5T_NATIVE_DOUBLE, diagsGroup.dataspace_id,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        diagsGroup.dataset_id.push_back(dataset_id);
        H5Dclose( diagsGroup.dataset_id[iDiag] );
        H5Sclose( diagsGroup.dataspace_id );



        // Create dataset for angleDist
        diagsGroup.dims_global[2] = 90;
        diagsGroup.dims_global[1] = 2;
        diagsGroup.dims_global[0] = diag1D->n_species;

        iDiag++;
        diagsGroup.dataspace_id = H5Screate_simple(data_dims, diagsGroup.dims_global, NULL);
        h5_name = diagsGroup.dataset_stringName[iDiag].c_str();
        dataset_id = H5Dcreate2(diagsGroup.group_id, h5_name, H5T_NATIVE_DOUBLE, diagsGroup.dataspace_id,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        diagsGroup.dataset_id.push_back(dataset_id);
        H5Dclose( diagsGroup.dataset_id[iDiag] );
        H5Sclose( diagsGroup.dataspace_id );


        // Create dataset for particleNumber
        diagsGroup.dims_global[2] = diag1D->n_species;
        diagsGroup.dims_global[1] = 1;
        diagsGroup.dims_global[0] = 1;

        iDiag++;
        diagsGroup.dataspace_id = H5Screate_simple(data_dims, diagsGroup.dims_global, NULL);
        h5_name = diagsGroup.dataset_stringName[iDiag].c_str();
        dataset_id = H5Dcreate2(diagsGroup.group_id, h5_name, H5T_NATIVE_INT, diagsGroup.dataspace_id,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        diagsGroup.dataset_id.push_back(dataset_id);
        H5Dclose( diagsGroup.dataset_id[iDiag] );
        H5Sclose( diagsGroup.dataspace_id );


        // Create dataset for kineticEnergy
        diagsGroup.dims_global[2] = diag1D->n_species;
        diagsGroup.dims_global[1] = 1;
        diagsGroup.dims_global[0] = 1;

        iDiag++;
        diagsGroup.dataspace_id = H5Screate_simple(data_dims, diagsGroup.dims_global, NULL);
        h5_name = diagsGroup.dataset_stringName[iDiag].c_str();
        dataset_id = H5Dcreate2(diagsGroup.group_id, h5_name, H5T_NATIVE_DOUBLE, diagsGroup.dataspace_id,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        diagsGroup.dataset_id.push_back(dataset_id);
        H5Dclose( diagsGroup.dataset_id[iDiag] );
        H5Sclose( diagsGroup.dataspace_id );

        diagsGroup.status = H5Gclose( diagsGroup.group_id );
    }
    diagsGroup.dataset_id.resize( diagsGroup.dataset_stringName.size() );
}


void SmileiIO_Cart1D::createPartsPattern( PicParams& params, ElectroMagn* fields, vector<Species*>& vecSpecies )
{

    // For particles, size ofdims_global should be 5: dims_global[nx][ny][nz][nvelocity][ntime]
    // But to be simple, the size is set 4, nz dimension is deleted.
    ptclsGroup.dims_global[2] = vx_dim;
    ptclsGroup.dims_global[1] = params.n_space_global[0];
    ptclsGroup.dims_global[0] = 1;

    ptclsGroup.ndims_[0] = ptclsGroup.dims_global[0];
    ptclsGroup.ndims_[1] = ptclsGroup.dims_global[1];
    ptclsGroup.ndims_[2] = ptclsGroup.dims_global[2];

    ptclsGroup.offset[0] = 0;
    ptclsGroup.offset[1] = 0;
    ptclsGroup.offset[2] = 0;

    ptclsGroup.stride[0] = 1;
    ptclsGroup.stride[1] = 1;
    ptclsGroup.stride[2] = 1;

    ptclsGroup.block[0] = 1;
    ptclsGroup.block[1] = 1;
    ptclsGroup.block[2] = 1;

    ptclsGroup.aDims = 3;

    createPartsGroup(vecSpecies);


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
    const char* h5_name;
    int data_dims = 3;
    Diagnostic1D* diag1D = static_cast<Diagnostic1D*>(diag);
    if(params.is_calVDF)
    {
        calVDF( params, smpi, fields, vecSpecies, itime);
    }
    if( itime % params.dump_step == 0 && smpi->isMaster() )
    {
        ndims_t = itime / params.dump_step - 1;
        long long ndims_t_temp = ndims_t;

        // create file at current output step
        data_file_name = "data/data" + to_string(ndims_t_temp) + ".h5";
        data_file_id = H5Fcreate( data_file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        // =============write fields============================================
        fieldsGroup.group_id = H5Gcreate(data_file_id, "/Fields", H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);
        for(int i = 0; i < fieldsGroup.dataset_stringName.size(); i++)
        {
            fieldsGroup.dataspace_id = H5Screate_simple(data_dims, fieldsGroup.dims_global, NULL);
            h5_name = fieldsGroup.dataset_stringName[i].c_str();
            fieldsGroup.dataset_id[i] = H5Dcreate2(fieldsGroup.group_id, h5_name, H5T_NATIVE_DOUBLE, fieldsGroup.dataspace_id,
                                                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            fieldsGroup.status = H5Dwrite(fieldsGroup.dataset_id[i], H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, fieldsGroup.dataset_data[i]);
            fieldsGroup.status = H5Sclose(fieldsGroup.dataspace_id);
            fieldsGroup.status = H5Dclose(fieldsGroup.dataset_id[i]);
        }
        fieldsGroup.status = H5Gclose(fieldsGroup.group_id);


/*
        // ==============write particle velocity distribution function=========================
        ptclsGroup.group_id = H5Gopen(data_file_id, "/VDF", H5P_DEFAULT);
        for(int i = 0; i < ptclsGroup.dataset_stringName.size(); i++)
        {
            h5_name = ptclsGroup.dataset_stringName[i].c_str();
            ptclsGroup.dataset_id[i] = H5Dopen( ptclsGroup.group_id, h5_name, H5P_DEFAULT);
        }

        ptclsGroup.count[0]  = ptclsGroup.dims_global[0];
        ptclsGroup.count[1]  = ptclsGroup.dims_global[1];
        ptclsGroup.count[2]  = ptclsGroup.dims_global[2];
        for(int i = 0; i < ptclsGroup.dataset_id.size(); i++)
        {
            ptclsGroup.memspace_id = H5Screate_simple (data_dims, ptclsGroup.count, NULL);
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
*/

        // ==============write Diagnostic: particleFlux, heatFlux and angleDist============
        // ==============particleNumber, kineticEnergy                         ============
        createDiagsPattern(params, diag1D);
        diagsGroup.group_id = H5Gopen(data_file_id, "/Diagnostic", H5P_DEFAULT);
        for(int i = 0; i < diagsGroup.dataset_stringName.size(); i++)
        {
            h5_name = diagsGroup.dataset_stringName[i].c_str();
            diagsGroup.dataset_id[i] = H5Dopen( diagsGroup.group_id, h5_name, H5P_DEFAULT);
        }

        for(int ispec = 0; ispec < diag1D->n_species; ispec++)
        {
            for(int iDirection = 0; iDirection < 2; iDirection++)
            {
                diagsGroup.offset[0] = 0;
                diagsGroup.offset[1] = ispec;
                diagsGroup.offset[2] = iDirection;
                diagsGroup.count[0]  = 1;
                diagsGroup.count[1]  = 1;
                diagsGroup.count[2]  = 1;

                // particleFlux
                diagsGroup.memspace_id = H5Screate_simple (data_dims, diagsGroup.count, NULL);
                diagsGroup.dataspace_id = H5Dget_space (diagsGroup.dataset_id[0]);
                diagsGroup.status = H5Sselect_hyperslab (diagsGroup.dataspace_id, H5S_SELECT_SET, diagsGroup.offset,
                                                  diagsGroup.stride, diagsGroup.count, diagsGroup.block);
                diagsGroup.status = H5Dwrite (diagsGroup.dataset_id[0], H5T_NATIVE_DOUBLE, diagsGroup.memspace_id,
                                     diagsGroup.dataspace_id, H5P_DEFAULT, &(diag1D->particleFlux[ispec][iDirection]) );

                diagsGroup.status = H5Sclose (diagsGroup.memspace_id);
                diagsGroup.status = H5Sclose (diagsGroup.dataspace_id);

                // heatFlux
                diagsGroup.memspace_id = H5Screate_simple (data_dims, diagsGroup.count, NULL);
                diagsGroup.dataspace_id = H5Dget_space (diagsGroup.dataset_id[1]);
                diagsGroup.status = H5Sselect_hyperslab (diagsGroup.dataspace_id, H5S_SELECT_SET, diagsGroup.offset,
                                                  diagsGroup.stride, diagsGroup.count, diagsGroup.block);
                diagsGroup.status = H5Dwrite (diagsGroup.dataset_id[1], H5T_NATIVE_DOUBLE, diagsGroup.memspace_id,
                                     diagsGroup.dataspace_id, H5P_DEFAULT, &(diag1D->heatFlux[ispec][iDirection]) );

                diagsGroup.status = H5Sclose (diagsGroup.memspace_id);
                diagsGroup.status = H5Sclose (diagsGroup.dataspace_id);


                diagsGroup.offset[0] = ispec;
                diagsGroup.offset[1] = iDirection;
                diagsGroup.offset[2] = 0;
  
                diagsGroup.count[0]  = 1;
                diagsGroup.count[1]  = 1;
                diagsGroup.count[2]  = 90;

                // angleDist
                diagsGroup.memspace_id = H5Screate_simple (data_dims, diagsGroup.count, NULL);
                diagsGroup.dataspace_id = H5Dget_space (diagsGroup.dataset_id[2]);
                diagsGroup.status = H5Sselect_hyperslab (diagsGroup.dataspace_id, H5S_SELECT_SET, diagsGroup.offset,
                                                  diagsGroup.stride, diagsGroup.count, diagsGroup.block);
                diagsGroup.status = H5Dwrite (diagsGroup.dataset_id[2], H5T_NATIVE_DOUBLE, diagsGroup.memspace_id,
                                     diagsGroup.dataspace_id, H5P_DEFAULT, &(diag1D->angleDist[ispec][iDirection][0]) );

                diagsGroup.status = H5Sclose(diagsGroup.memspace_id);
                diagsGroup.status = H5Sclose(diagsGroup.dataspace_id);

            }
        }

        diagsGroup.offset[0] = 0;
        diagsGroup.offset[1] = 0;
        diagsGroup.offset[2] = 0;

        diagsGroup.count[0]  = 1;
        diagsGroup.count[1]  = 1;
        diagsGroup.count[2]  = diag1D->n_species;

        // particleNumber
        diagsGroup.memspace_id = H5Screate_simple (data_dims, diagsGroup.count, NULL);
        diagsGroup.dataspace_id = H5Dget_space (diagsGroup.dataset_id[3]);
        diagsGroup.status = H5Sselect_hyperslab (diagsGroup.dataspace_id, H5S_SELECT_SET, diagsGroup.offset,
                                          diagsGroup.stride, diagsGroup.count, diagsGroup.block);
        diagsGroup.status = H5Dwrite (diagsGroup.dataset_id[3], H5T_NATIVE_INT, diagsGroup.memspace_id,
                             diagsGroup.dataspace_id, H5P_DEFAULT, &(diag1D->particleNumber[0]) );

        diagsGroup.status = H5Sclose (diagsGroup.memspace_id);
        diagsGroup.status = H5Sclose (diagsGroup.dataspace_id);


        // kineticEnergy
        diagsGroup.memspace_id = H5Screate_simple (data_dims, diagsGroup.count, NULL);
        diagsGroup.dataspace_id = H5Dget_space (diagsGroup.dataset_id[4]);
        diagsGroup.status = H5Sselect_hyperslab (diagsGroup.dataspace_id, H5S_SELECT_SET, diagsGroup.offset,
                                          diagsGroup.stride, diagsGroup.count, diagsGroup.block);
        diagsGroup.status = H5Dwrite (diagsGroup.dataset_id[4], H5T_NATIVE_DOUBLE, diagsGroup.memspace_id,
                             diagsGroup.dataspace_id, H5P_DEFAULT, &(diag1D->kineticEnergy[0]) );

        diagsGroup.status = H5Sclose (diagsGroup.memspace_id);
        diagsGroup.status = H5Sclose (diagsGroup.dataspace_id);
        


        // close dataset, group, and so on
        for(int i = 0; i < diagsGroup.dataset_stringName.size(); i++)
        {
            diagsGroup.status = H5Dclose(diagsGroup.dataset_id[i]);
        }
        diagsGroup.status = H5Gclose( diagsGroup.group_id );
        diagsGroup.dataset_stringName.clear();
        diagsGroup.dataset_id.clear();

        // close hdf5 file
        status = H5Fclose(data_file_id);
    }


} // END write
