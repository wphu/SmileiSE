/*
 * SmileIO_Cart2D.h
 *
 *  Created on: 3 juil. 2013
 */
#ifndef SMILEIO_CART2D_H
#define SMILEIO_CART2D_H

#include <string>
#include <vector>

#include "SmileiIO.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class SmileiIO_Cart2D
//  --------------------------------------------------------------------------------------------------------------------
class SmileiIO_Cart2D : public SmileiIO {
public:
    //! Create // HDF5 environment
    SmileiIO_Cart2D( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies );
    //! Destructor for SmileiIO
    ~SmileiIO_Cart2D();

    //! Build memory and file space for // HDF5 write/read
    void createPattern( PicParams& params, SmileiMPI* smpi, ElectroMagn* fields );

private:



};

#endif /* SMILEIO_CART2D_H_ */
