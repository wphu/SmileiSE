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

SmileiIO_Cart1D::SmileiIO_Cart1D( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies  )
: SmileiIO( params, smpi )
{
    initVDF(params, smpi, fields, vecSpecies);
    if(smpi->isMaster()) {
        createFieldsPattern(params, smpi, fields);
        createPartsPattern(params, smpi, fields, vecSpecies);
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


    fieldsGroup.group_id = H5Gcreate(global_file_id_, "/1d_global", H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);
    H5Gcreate(global_file_id_, "/2d_global", H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);

    /* Create a datagroup attribute. */
    dims = 4;
    //> dataspace is to descript the structure of data: the number of data dimension and the size of each dimension
    //> the first parameter rank=1: is the number of dimensions used in the dataspace
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
    for(int i = 0; i < fields->rho_s.size(); i++)
    {
        addField(fields->rho_s_global[i]);
        addField(fields->rho_s_global_avg[i]);

    }

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
    /* Close the attribute. */
    ptclsGroup.status = H5Aclose(ptclsGroup.attribute_id);



    Species *s;
    Particles *p;


    for(int isp=0; isp<vx_VDF.size(); isp++)
    {
        s = vecSpecies[isp];
        p = &(s->particles);

        string fieldName = "VDF_" + s->species_param.species_type;
        const char* name = fieldName.c_str();
        ptclsGroup.dataset_name.push_back(name);

        /* Create the data space for the dataset. */
        ptclsGroup.dataspace_id = H5Screate_simple(4, ptclsGroup.dims_global, NULL);
        int dataset_size = ptclsGroup.dataset_id.size();
        hid_t id = H5Dcreate2(ptclsGroup.group_id, name, H5T_NATIVE_DOUBLE, ptclsGroup.dataspace_id,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        ptclsGroup.dataset_id.push_back(id);

        ptclsGroup.dataset_data.push_back(vx_VDF_global[isp]->data_);

        // for VDF_tot_ the dimension is same with VDF_, but only [nx=0] is meaningful
        fieldName = "VDF_tot_" + s->species_param.species_type;
        name = fieldName.c_str();
        ptclsGroup.dataset_name.push_back(name);

        /* Create the data space for the dataset. */
        ptclsGroup.dataspace_id = H5Screate_simple(4, ptclsGroup.dims_global, NULL);
        dataset_size = ptclsGroup.dataset_id.size();
        id = H5Dcreate2(ptclsGroup.group_id, name, H5T_NATIVE_DOUBLE, ptclsGroup.dataspace_id,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        ptclsGroup.dataset_id.push_back(id);

        ptclsGroup.dataset_data.push_back(vx_VDF_tot_global[isp]->data_);


    }

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
