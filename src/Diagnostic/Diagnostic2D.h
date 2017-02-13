#ifndef DIAGNOSTIC2D_H
#define DIAGNOSTIC2D_H

#include "Diagnostic.h"
#include "PicParams.h"
#include "SmileiMPI.h"


class Diagnostic2D : public Diagnostic {

public :

    Diagnostic2D(PicParams& params, SmileiMPI* smpi);
    virtual ~Diagnostic2D() {};

    //! Runs the diag for all patches for local diags.
    virtual void run( SmileiMPI* smpi, vector<Species*>& vecSpecies, ElectroMagn* EMfields, int itime ) {};


protected :


};

#endif
