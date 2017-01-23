#ifndef DIAGNOSTIC2D_H
#define DIAGNOSTIC2D_H

#include "Diagnostic.h"
#include "PicParams.h"
#include "SmileiMPI.h"


class Diagnostic2D : public Diagnostic {

public :

    Diagnostic2D() {};
    virtual ~Diagnostic2D() {};


    //! Misc init.
    virtual void init(PicParams& params, SmileiMPI* smpi) {};

    //! Runs the diag for a given patch for global diags.
    virtual void run( int timestep ) {};

    //! Runs the diag for all patches for local diags.
    virtual void run( SmileiMPI* smpi, int timestep ) {};


protected :


};

#endif
