#ifndef GRID_H
#define GRID_H

#include <cmath>

#include <vector>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>

#include "Tools.h"
#include "PicParams.h"


//! Class Grid: generic class allowing to define complex geometry and boundary, also used
//! to decide if the particle hit the wall
class Grid
{

public:



    //! Constructor for Grid: with no input argument
    Grid(){};

    Grid(PicParams &params){};

    //! Destructor for Grid
    virtual ~Grid() {
        ;
    } ;

    //! Virtual method used to allocate Grid
    virtual void allocateDims(){};
    virtual void geometry(){};
    virtual void computeNcp(){};

    //! vector containing the dimensions of the Grid
    //! \todo private/friend/modify SmileiMPI* (JD)
    std::vector<int> dims_;
    std::vector<int> globalDims_;

    //! returns the dimension of the Grid
    inline std::vector<int> dims () {return dims_;}
    //! All arrays may be viewed as a 1D array
    //! Linearized diags


    //! pointer to the linearized array
    int* iswall_;

    int* iswall_global_;
    int* bndr_global_;
    double* bndrVal_global_;
    //! The number of the current point in the discrete Poisson Eqution left coefficient matrix
    int* numcp_global_;
    //>>>total number of numcp_global_ points
    int ncp;
    std::vector<int> dims_source;
    int nx,ny,nz;


protected:

private:

};

#endif
