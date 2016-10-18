/*
 * SmileIO_Cart2D.h
 *
 *  Created on: 3 juil. 2013
 */
#ifndef SMILEIO_CART1D_H
#define SMILEIO_CART1D_H

#include <string>
#include <vector>

#include "SmileiIO.h"
#include "Array4D.h"

using namespace std;

//  --------------------------------------------------------------------------------------------------------------------
//! Class SmileiIO_Cart1D
//  --------------------------------------------------------------------------------------------------------------------
class SmileiIO_Cart1D : public SmileiIO {
public:
    //! Create // HDF5 environment
    SmileiIO_Cart1D( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies );
    //! Destructor for SmileiIO
    ~SmileiIO_Cart1D();

    //! Build memory and file space for // HDF5 write/read
    void createFieldsPattern( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields );

    // Create particles h5 file pattern
    void createPartsPattern( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies );

    void initVDF( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies );
    // calculate velocity distribution function
    void calVDF( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies);

    vector<Array4D*> vx_VDF;
    vector<Array4D*> vx_VDF_global;
    double vxMin,vxMax;
    double vx_d;
    int vx_dim;

private:


};

#endif /* SMILEIO_CART1D_H_ */
