#include "Grid2D.h"

#include <iostream>
#include <vector>
#include <cstring>
#include <iomanip>
#include <fstream>

using namespace std;



// ---------------------------------------------------------------------------------------------------------------------
// Creators for Grid2D
// ---------------------------------------------------------------------------------------------------------------------

// with no input argument
Grid2D::Grid2D() : Grid()
{


}

// with the dimensions as input argument
Grid2D::Grid2D(
    PicParams &params,
    string grid_type,
    string gap_kind,
    int ny_source_temp,
    int ny_gapHeight_temp,
    int nx_gapWeight_temp,
    double potential_wall_temp):
    Grid(params)
{
    gridType = grid_type;
    gapKind = gap_kind;
    ny_source = ny_source_temp;
    ny_gapHeight = ny_gapHeight_temp;
    nx_gapWeight = nx_gapWeight_temp;
    potential_wall = potential_wall_temp;

    // number of nodes of the grid in the x-direction
    dims_.resize(2);
    globalDims_.resize(2);

    dims_[0] = params.n_space[0]+1+2*params.oversize[0];
    dims_[1] = params.n_space[1]+1+2*params.oversize[1];

    globalDims_[0]=params.n_space_global[0]+1;
    globalDims_[1]=params.n_space_global[1]+1;

    nx=globalDims_[0];
    ny=globalDims_[1];

    allocateDims();
    if(gridType == "rectangle")
    {
            geometry();
    }
    else if(gridType == "gap")
    {
        geometry_gap();
    }
    computeNcp();
}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Grid2D
// ---------------------------------------------------------------------------------------------------------------------
Grid2D::~Grid2D()
{

}


void Grid2D::allocateDims( )
{
    iswall_             = new int[dims_[0]*dims_[1]];
    iswall_global_      = new int[globalDims_[0]*globalDims_[1]];
    bndr_global_        = new int[globalDims_[0]*globalDims_[1]];
    bndrVal_global_     = new double[globalDims_[0]*globalDims_[1]];
    numcp_global_       = new int[globalDims_[0]*globalDims_[1]];

    iswall_2D           =new int*[dims_[0]];
    for (int i=0; i<dims_[0]; i++)  {
        iswall_2D[i] = iswall_ + i*dims_[1];
        for (int j=0; j<dims_[1]; j++) iswall_2D[i][j] = 0;
    }

    iswall_global_2D    =new int*[globalDims_[0]];
    for (int i=0; i<globalDims_[0]; i++)  {
        iswall_global_2D[i] = iswall_global_ + i*globalDims_[1];
        for (int j=0; j<globalDims_[1]; j++) iswall_global_2D[i][j] = 0;
    }

    bndr_global_2D      =new int*[globalDims_[0]];
    for (int i=0; i<globalDims_[0]; i++)  {
        bndr_global_2D[i] = bndr_global_ + i*globalDims_[1];
        for (int j=0; j<globalDims_[1]; j++) bndr_global_2D[i][j] = 0;
    }

    bndrVal_global_2D   =new double*[globalDims_[0]];
    for (int i=0; i<globalDims_[0]; i++)  {
        bndrVal_global_2D[i] = bndrVal_global_ + i*globalDims_[1];
        for (int j=0; j<globalDims_[1]; j++) bndrVal_global_2D[i][j] = 0.0;
    }

    numcp_global_2D    =new int*[globalDims_[0]];
    for (int i=0; i<globalDims_[0]; i++)  {
        numcp_global_2D[i] = numcp_global_ + i*globalDims_[1];
        for (int j=0; j<globalDims_[1]; j++) numcp_global_2D[i][j] = 0;
    }

}


//>>>no gap geometry, with source in x direction
void Grid2D::geometry( )
{


    dims_source.resize(2);
    dims_source[0]=0;
    dims_source[1]=0;

    // iswall_global_2D is for particle moving, if four points of one grid is 1, then the grid is wall
    for(int i=0; i<nx; i++)
      for(int j=0; j<ny; j++){
        iswall_global_2D[i][j]=0;
      }

    for(int i=0; i<nx; i++){
      iswall_global_2D[i][0]=1;
      iswall_global_2D[i][ny-1]=1;
    }

    for(int j=0; j<ny; j++){
      iswall_global_2D[0][j]=1;
      iswall_global_2D[nx-1][j]=1;
    }

    //>>>struct boundary condition
    // bndr* is for electric potential solving
    for(int i=0; i<nx; i++)
      for(int j=0; j<ny; j++){
        bndr_global_2D[i][j]=0;
      }

    // 5 is the particle source region, usually the region has no need to solve the electric potential
    // The electric potential in the region can be set to some constant
    for(int i=0; i<dims_source[0]; i++)
      for(int j=0; j<ny; j++){
        bndr_global_2D[i][j]=5;
      }

    for(int i=0; i<nx; i++)
      for(int j=0; j<dims_source[1]; j++){
        bndr_global_2D[i][j]=5;
      }

    /*
    // 8 is the periodic boundary condition
    for(int i=dims_source[0]; i<nx; i++){
      bndr_global_2D[i][0]=8;
      bndr_global_2D[i][ny-1]=8;
    }
    */

    for(int i=dims_source[0]; i<nx; i++){
      bndr_global_2D[i][0] = 1;
      bndrVal_global_2D[i][0] = 0.0;

      bndr_global_2D[i][ny-1] = 1;
      bndrVal_global_2D[i][ny-1] = 0.0;
    }


    // 1 is the Dirchlet boundary condition
    for(int j=0; j<ny; j++){
      bndr_global_2D[dims_source[0]][j]=1;
      bndrVal_global_2D[dims_source[0]][j]=0.0;

      bndr_global_2D[nx-1][j]=1;
      bndrVal_global_2D[nx-1][j]=0.0;
    }


}


//>>>classical gap geometry, with source in x direction
void Grid2D::geometry_gap( )
{

    ofstream isWall;
    ofstream bndr;
    ofstream bndrVal;

    // electric potential at the wall surface
    double val1 = potential_wall;

    // electric potential at source region
    double val1_source = 0.0;


    // iswall_global_2D is for particle moving, if four points of one grid is 1, then the grid is wall
    for(int i=0; i<nx; i++)
      for(int j=0; j<ny; j++){
        iswall_global_2D[i][j]=0;
      }

    for(int i=0; i<nx; i++){
      iswall_global_2D[i][0]=1;
      iswall_global_2D[i][ny-1]=1;
    }

    for(int j=0; j<ny; j++){
      iswall_global_2D[0][j]=1;
      iswall_global_2D[nx-1][j]=1;
    }

    for(int i = 0; i < 0.5*nx - 0.5*nx_gapWeight; i++)
    {
        for(int j = 0; j < ny_gapHeight; j++)
        {
            iswall_global_2D[i][j]=1;
        }
    }

    for(int i =  0.5*nx + 0.5*nx_gapWeight; i < nx; i++)
    {
        for(int j = 0; j < ny_gapHeight; j++)
        {
            iswall_global_2D[i][j]=1;
        }
    }


    //==============struct boundary condition===================================
    // bndr* is for electric potential solving
    for(int i=0; i<nx; i++)
      for(int j=0; j<ny; j++){
        bndr_global_2D[i][j]=0;
      }

    // 5 is the particle source region or wall region, usually the region has no need to solve the electric potential
    // The electric potential in the region can be set to some constant
    for(int i=0; i<nx; i++)
    {
        for(int j = ny - ny_source; j < ny; j++)
        {
            bndr_global_2D[i][j] = 5;
            bndrVal_global_2D[i][j] = val1_source;
        }
    }

    // 8 is the periodic boundary condition
    for(int j = 0; j < ny - ny_source; j++){
      bndr_global_2D[0][j]=8;
      bndr_global_2D[nx-1][j]=8;
    }

    // define the source region surface boudary
    for(int i = 0; i < nx; i++){
      bndr_global_2D    [i][ny - ny_source -1] = 1;
      bndrVal_global_2D [i][ny - ny_source -1] = val1_source;
    }

    // define the wall surface boundary
    for(int j = 0; j < ny_gapHeight; j++)
    {
        for(int i = 0; i < nx; i++)
        {
            if( (i == 0 || i == nx - 1) && j == ny_gapHeight -1 ) {
                bndr_global_2D[i][j] = 1;
                bndrVal_global_2D[i][j] = val1;
            }
            else if( (i == 0 || i == nx - 1) && j < ny_gapHeight -1 ) {
                bndr_global_2D[i][j] = 5;
                bndrVal_global_2D[i][j] = val1;
            }
            else if( j == 0 && ( iswall_global_2D[i][j+1] == 0 ||
            (i != 0 && iswall_global_2D[i-1][j+1] == 0) || (i != nx-1 && iswall_global_2D[i+1][j+1] == 0) ) ) {
                bndr_global_2D[i][j] = 1;
                bndrVal_global_2D[i][j] = val1;
            }
            else if( j == 0 && iswall_global_2D[i][j+1] == 1 ) {
                bndr_global_2D[i][j] = 5;
                bndrVal_global_2D[i][j] = val1;
            }
            else if( iswall_global_2D[i][j] == 1 && ( iswall_global_2D[i-1][j] == 0 || iswall_global_2D[i][j-1] == 0
             || iswall_global_2D[i+1][j] == 0 || iswall_global_2D[i][j+1] == 0 ) ) {
                bndr_global_2D[i][j] = 1;
                bndrVal_global_2D[i][j] = val1;
            }
            else if( iswall_global_2D[i][j] == 1 ) {
                bndr_global_2D[i][j] = 5;
                bndrVal_global_2D[i][j] = val1;
            }
        }
    }


    // output data grid data
    isWall.open("isWall.txt");
    bndr.open("bndr.txt");
    bndrVal.open("bndrVal.txt");
    for(int i=0; i<nx; i++)
    {
        for(int j = 0; j < ny; j++)
        {
            isWall<<iswall_global_2D[i][j];
            bndr<<bndr_global_2D[i][j];
            bndrVal<<bndrVal_global_2D[i][j];
        }
        isWall<<endl;
        bndr<<endl;
        bndrVal<<endl;
    }
    isWall.close();
    bndr.close();
    bndrVal.close();

}



void Grid2D::computeNcp(){

    ncp=0;
    for(int i=0; i<nx; i++)
      for(int j=0; j<ny; j++){
        if( bndr_global_2D[i][j]==0 || bndr_global_2D[i][j]==1
        || bndr_global_2D[i][j]==2 || bndr_global_2D[i][j]==8){
          ncp++;
          numcp_global_2D[i][j]=ncp-1;

        }
      }


}
