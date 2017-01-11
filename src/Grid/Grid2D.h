#ifndef GRID2D_H
#define GRID2D_H

#include <cmath>

#include <vector>
#include <string>

#include "Grid.h"

using namespace std;

//! class Grid2D used to defined a 2d vector
class Grid2D : public Grid
{

public:
    //! Constructor for Grid2D: no input argument
    Grid2D();

    //! Constructor for Grid2D: with the vector dimension as input argument
    Grid2D(
        PicParams &params,
        string grid_type,
        string gap_kind,
        int ny_source_temp,
        int ny_gapHeight_temp,
        int nx_gapWeight_temp,
        double potential_wall_temp);

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



    // Tomakak divertor gap geometry Parameters
    std::string gapKind;
    int ny_source;
    int ny_gapHeight;
    int nx_gapWeight;
    double potential_wall;


private:
    //!\todo{Comment what are these stuffs (MG for JD)}
    //double *data_2D;
};

#endif
