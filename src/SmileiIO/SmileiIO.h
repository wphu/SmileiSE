/*
 * SmileiIO.h
 *
 *  Created on: 3 juil. 2013
 */

#ifndef SMILEIIO_H
#define SMILEIIO_H

#include <string>
#include <vector>

#include <Field.h>
#include <hdf5.h>
#include <Tools.h>
#include "Array4D.h"

class PicParams;
class InputData;
class SmileiMPI;
class ElectroMagn;
class Field;
class Species;
class Diagnostic;

#include <csignal>

using namespace std;

//  --------------------------------------------------------------------------------------------------------------------
//! Class SmileiIO
//  --------------------------------------------------------------------------------------------------------------------
class SmileiIO {
public:
    //! Create // HDF5 environment
    //! @see global_file_id_
    //! @see global_file_id_avg
    SmileiIO( PicParams& params, SmileiMPI* smpi );
    //! Destructor for SmileiIO
    virtual ~SmileiIO();

    void addField(Field* field);
    //! Basic write field on its own file (debug)
    virtual void write( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies, Diagnostic* diag){};
    virtual void calVDF( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies){};


    //> Id of "Fields_global.h5", contains global fields, such as potential, rho ...
    hid_t       global_file_id_;
    herr_t      status;

    double* data_;

    //! Space dimension of a particle
    unsigned int nDim_particle;

    //> dimensions of time, which means total timestep number to output
    int ndims_t;

    // H5Group is for fields and particles( velocity distribution funtion )
    struct H5Group
    {
        hid_t       group_id;
        hid_t       attribute_id;
        herr_t      status;

        // dataspace_id is temporary to create dataset_id, corresponse to the whole dataset structure
        hid_t       dataspace_id;
        // memspace_id is temporary to write subset of a dataset, descript the structure of the subset
        hid_t       memspace_id;
        // ndims_ is used to output the dimensitons as attribute
        int ndims_[4];

        std::vector<hid_t> dataset_id;
        std::vector<std::string> dataset_stringName;
        std::vector<const char*> dataset_name;
        std::vector<double*> dataset_data;

        //> data dimensions to be outputed: t, z, y, x
        //> [1] = 1 for 2d; [1] = [2] = 1 for 1d
        // For different dataset, the dims_global may be different
        hsize_t     dims_global[4];

        //> parameters for outputing fields in one timestep, used by H5Sselect_hyperslab hdf5 method
        /*
        The offset or start array specifies the offset of the starting element of the specified hyperslab.
        The count array determines how many blocks to select from the dataspace in each dimension. If the block size for a dimension is one then the count is the number of elements along that dimension.
        The stride array allows you to sample elements along a dimension. For example, a stride of one (or NULL) will select every element along a dimension, a stride of two will select every other element, and a stride of three will select an element after every two elements.
        The block array determines the size of the element block selected from a dataspace. If the block size is one or NULL then the block size is a single element in that dimension
        */
        hsize_t     count[4];              /* size of subset in the file */
        hsize_t     offset[4];             /* subset offset in the file */
        hsize_t     stride[4];
        hsize_t     block[4];
    };

    H5Group fieldsGroup;
    H5Group ptclsGroup;     //H5 particles Group
    H5Group diagsGroup;     //Diagnostic Group


    std::vector<Array4D*> vx_VDF;
    std::vector<Array4D*> vx_VDF_global;
    std::vector<Array4D*> vx_VDF_tot_global;     // [1][1][1][nv] velocity distribution function for particles in all cells
    double vxMin,vxMax;
    double vx_d;
    int vx_dim;



private:
    //! incremental number of times we've done a dump
    unsigned int dump_times;

    //! name of the fields to dump
    std::vector<std::string> fieldsToDump;





};

#endif /* SMILEI_OUTPUT_H_ */
