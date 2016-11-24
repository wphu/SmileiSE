#ifndef GRIDFACTORY_H
#define GRIDFACTORY_H

#include "Grid.h"
#include "Grid2D.h"

#include "PicParams.h"

#include "Tools.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class SmileiMPIFactory
//  --------------------------------------------------------------------------------------------------------------------
class GridFactory {
public:
    //  --------------------------------------------------------------------------------------------------------------------
    //! Create appropriate MPI environment for the geometry
    //! \param params : Parameters
    //! \param smpiData : Initial MPI environment (data broadcast)
    //  --------------------------------------------------------------------------------------------------------------------
    static Grid* create(PicParams& params) {
        Grid* grid = NULL;
        MESSAGE(1, "Geometry:" << params.geometry);
        if ( params.geometry == "1d3v" ) {

        }
        else if ( params.geometry == "2d3v" ) {
            grid = new Grid2D(params);
        }
        else {
            ERROR( "Geometry " << params.geometry << " not implemented" );
        }

        return grid;
    }

};

#endif
