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

class PicParams;
class InputData;
class SmileiMPI;
class ElectroMagn;
class Field;
class Species;

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
    void write( PicParams& params, ElectroMagn* fields, SmileiMPI* smpi, vector<Species*>& vecSpecies);
    virtual void calVDF( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies){};


    //> Id of "Fields_global.h5", contains global fields, such as potential, rho ...
    hid_t global_file_id_;

    double* data_;

    //! Space dimension of a particle
    unsigned int nDim_particle;

    //> dimensions of time, which means total timestep number to output
    int ndims_t;

    struct H5Group
    {
        hid_t       group_id;
        hid_t       dataspace_id;
        hid_t       attribute_id;
        hid_t       memspace_id;
        herr_t      status;

        std::vector<hid_t> dataset_id;
        std::vector<const char*> dataset_name;
        std::vector<double*> dataset_data;

        //> data dimensions to be outputed: t, z, y, x
        //> ndims_[1] = 1 for 2d; ndims_[1] = ndims_[2] = 1 for 1d
        int ndims_[4];
        hsize_t     dims_global[4];

        //> parameters for outputing fields in one timestep, used by H5Sselect_hyperslab hdf5 method
        hsize_t     count[4];              /* size of subset in the file */
        hsize_t     offset[4];             /* subset offset in the file */
        hsize_t     stride[4];
        hsize_t     block[4];
    };

    H5Group fieldsGroup;
    H5Group ptclsGroup;     //H5 particles Group


private:
    //! incremental number of times we've done a dump
    unsigned int dump_times;

    //! name of the fields to dump
    std::vector<std::string> fieldsToDump;





};

#endif /* SMILEI_OUTPUT_H_ */
