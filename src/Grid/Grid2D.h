#ifndef GRID2D_H
#define GRID2D_H

#include <cmath>

#include <vector>

#include "Grid.h"

//! class Grid2D used to defined a 2d vector
class Grid2D : public Grid
{

public:
    //! Constructor for Grid2D: no input argument
    Grid2D();

    //! Constructor for Grid2D: with the vector dimension as input argument
    Grid2D(PicParams &params);

    //! Destructor for Grid2D
    ~Grid2D();

    //! Method used to allocate a Grid2D
    void allocateDims();
    void geometry();
    void geometry_gap();
    void computeNcp();
    int **iswall_2D;
  	int **iswall_global_2D;
    int **bndr_global_2D;
    double **bndrVal_global_2D;
    //! The number of the current point in the discrete Poisson Eqution left coefficient matrix
    int **numcp_global_2D;

    //>>>temp variable



private:
    //!\todo{Comment what are these stuffs (MG for JD)}
    //double *data_2D;
};

#endif
