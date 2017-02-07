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

    static Diagnostic* create(PicParams& params, SmileiMPI* smpi) {
        Diagnostic* Diag;

        return Diag;
    } // END createGlobalDiagnostics


};

#endif
