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
    virtual void createFieldsGroup( ElectroMagn* fields );
    virtual void createPartsGroup( vector<Species*>& vecSpecies );


    //! Basic write field on its own file (debug)
    virtual void write( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies, Diagnostic* diag, int itime){};
    //virtual void calVDF( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies){};

    // store Particles to restart
    virtual void storeP(PicParams& params, SmileiMPI* smpi, vector<Species*>& vecSpecies, int itime );

    // reload Particles to restart
    virtual void reloadP(PicParams& params, SmileiMPI* smpi, vector<Species*>& vecSpecies );

    // store Particles to restart
    virtual void endStoreP(PicParams& params, SmileiMPI* smpi, vector<Species*>& vecSpecies, int itime );

    string data_file_name;

    // Id of "Fields_global.h5", contains global fields, such as potential, rho ...
    hid_t       data_file_id;
    herr_t      status;

    double* data_;
    int* data_int;

    // Space dimension of a particle
    unsigned int nDim_particle;

    // dimensions of time, which means total timestep number to output
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
        int ndims_[3];
        // aDims is for For attribute
        hsize_t     aDims;

        std::vector<hid_t>          dataset_id;
        std::vector<std::string>    dataset_stringName;
        std::vector<double*>        dataset_data;

        //> data dimensions to be outputed: t, z, y, x
        //> [1] = 1 for 2d; [1] = [2] = 1 for 1d
        // For different dataset, the dims_global may be different
        hsize_t     dims_global[3];

        //> parameters for outputing fields in one timestep, used by H5Sselect_hyperslab hdf5 method
        /*
        The offset or start array specifies the offset of the starting element of the specified hyperslab.
        The count array determines how many blocks to select from the dataspace in each dimension. If the block size for a dimension is one then the count is the number of elements along that dimension.
        The stride array allows you to sample elements along a dimension. For example, a stride of one (or NULL) will select every element along a dimension, a stride of two will select every other element, and a stride of three will select an element after every two elements.
        The block array determines the size of the element block selected from a dataspace. If the block size is one or NULL then the block size is a single element in that dimension
        */
        hsize_t     count[3];              /* size of subset in the file */
        hsize_t     offset[3];             /* subset offset in the file */
        hsize_t     stride[3];
        hsize_t     block[3];
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

    // ======================= Some paramter for restart =======================
    // Id of "Restore000.h5", contains paritcles information to restart
    hid_t       restore_file_id_;
    string      restore_file_dir_;

    // If is_restart == 1, it means no need to create file "Fields_global.h5"
    int is_restart;

    int stepStart;
    int nSpecies;


private:
    //! incremental number of times we've done a dump
    unsigned int dump_times;

    //! name of the fields to dump
    std::vector<std::string> fieldsToDump;





};

#endif /* SMILEI_OUTPUT_H_ */
