#ifndef SMILEIIOFACTORY_H
#define SMILEIIOFACTORY_H

#include "SmileiIO.h"
#include "SmileiIO_Cart1D.h"
#include "SmileiIO_Cart2D.h"

#include "PicParams.h"
#include "SmileiMPI.h"
#include "Diagnostic.h"

#include "Tools.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class SmileiIOFactory
//  --------------------------------------------------------------------------------------------------------------------
class SmileiIOFactory {
public:
    //  --------------------------------------------------------------------------------------------------------------------
    //! Create appropriate IO environment for the geometry
    //! \param params : Parameters
    //! \param diag : Diagnostics
    //! \param smpi : MPI environment
    //  --------------------------------------------------------------------------------------------------------------------
    static SmileiIO* create(PicParams& params, SmileiMPI* smpi, ElectroMagn* fields, vector<Species*>& vecSpecies, Diagnostic* diag ) {
        SmileiIO* sio = NULL;
        if ( params.geometry == "1d3v" ) {
            sio = new  SmileiIO_Cart1D(params, smpi, fields, vecSpecies, diag);
        }
        else if ( params.geometry == "2d3v" ) {
            sio = new  SmileiIO_Cart2D(params, smpi, fields, vecSpecies, diag);
        }
        else {
            ERROR( "Geometry " << params.geometry << " not implemented" );
        }

        return sio;
    }

};

#endif
