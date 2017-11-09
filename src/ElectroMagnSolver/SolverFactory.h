#ifndef SOLVERFACTORY_H
#define SOLVERFACTORY_H

#include "EF_Solver1D_TDMA.h"
#include "EF_Solver1D_TDMA_imp.h"
#include "EF_Solver1D_GeneralThomas.h"
#include "MF_Solver1D_Yee.h"
#include "MF_Solver2D_Yee.h"
#include "MF_Solver2D_Cowan.h"
#include "EF_Solver2D_SLU.h"
#include "InputData.h"
#include "PicParams.h"

#include "Tools.h"

class SolverFactory {
public:
    static Solver* create(PicParams& params, InputData &ifile, Grid* grid, SmileiMPI* smpi) {
        Solver* solver = NULL;
        int nx_source_left;
        if ( params.geometry == "1d3v" ) {
            nx_source_left = 0.0; // default
            ifile.extract("nx_source_left", nx_source_left);
            if(params.method == "explicit")
            {
                if(params.solver_type == "GeneralThomas")
                {
                  solver = new EF_Solver1D_GeneralThomas(params, smpi, nx_source_left);
                }
                else
                {
                  solver = new EF_Solver1D_TDMA(params, smpi, nx_source_left);
                }

            }
            else if(params.method == "implicit")
            {
                solver = new EF_Solver1D_TDMA_imp(params, smpi, nx_source_left);
            }

        }
        else if ( params.geometry == "2d3v" ) {
	    //if ()
            //solver = new MF_Solver2D_Yee(params);
            solver = new EF_Solver2D_SLU(params, grid, smpi);
	    //elseif()
	    //solver = new MF_Solver1D_Cowan(params);
        }
        else {}

        return solver;
    }

};


#endif
