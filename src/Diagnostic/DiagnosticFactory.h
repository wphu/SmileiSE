#ifndef DIAGNOSTICFACTORY_H
#define DIAGNOSTICFACTORY_H

#include "Diagnostic1D.h"
#include "Diagnostic2D.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Create appropriate IO environment for the geometry
//! \param params : Parameters
//! \param smpi : MPI environment
//  --------------------------------------------------------------------------------------------------------------------

class DiagnosticFactory {
public:

    static Diagnostic* create(PicParams& params, SmileiMPI* smpi, ElectroMagn* EMfields) {
        Diagnostic* diag = NULL;
        if ( params.geometry == "1d3v" ) {
            diag = new Diagnostic1D(params, smpi, EMfields);
        }
        else if ( params.geometry == "2d3v" ) {
            diag = new Diagnostic2D(params, smpi, EMfields);
        }
        else {
            ERROR( "Unknwon geometry : " << params.geometry );
        }

        return diag;
    } // END createGlobalDiagnostics


};

#endif
