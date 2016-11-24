#include "Grid2D.h"

#include <iostream>
#include <vector>
#include <cstring>
#include <iomanip>


using namespace std;



// ---------------------------------------------------------------------------------------------------------------------
// Creators for Grid2D
// ---------------------------------------------------------------------------------------------------------------------

// with no input argument
Grid2D::Grid2D() : Grid()
{


}

// with the dimensions as input argument
Grid2D::Grid2D(PicParams &params):
Grid(params)
{
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
    geometry();
    computeNcp();
}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Grid2D
// ---------------------------------------------------------------------------------------------------------------------
Grid2D::~Grid2D()
{

}


void Grid2D::allocateDims( ){
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
void Grid2D::geometry( ){


    dims_source.resize(2);
    dims_source[0]=0;
    dims_source[1]=0;

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
    for(int i=0; i<nx; i++)
      for(int j=0; j<ny; j++){
        bndr_global_2D[i][j]=0;
      }

    for(int i=0; i<dims_source[0]; i++)
      for(int j=0; j<ny; j++){
        bndr_global_2D[i][j]=5;
      }

    for(int i=0; i<nx; i++)
      for(int j=0; j<dims_source[1]; j++){
        bndr_global_2D[i][j]=5;
      }

    for(int i=dims_source[0]; i<nx; i++){
      bndr_global_2D[i][0]=8;
      bndr_global_2D[i][ny-1]=8;
    }

    for(int j=0; j<ny; j++){
      bndr_global_2D[dims_source[0]][j]=1;
      bndrVal_global_2D[dims_source[0]][j]=0.0;

      bndr_global_2D[nx-1][j]=1;
      bndrVal_global_2D[nx-1][j]=0.0;
    }


}




void Grid2D::computeNcp(){

    ncp=0;
    for(int i=0; i<nx; i++)
      for(int j=0; j<ny; j++){
        if((iswall_global_2D[i][j]==0 && bndr_global_2D[i][j]!=5) || bndr_global_2D[i][j]==1
        || bndr_global_2D[i][j]==2 || bndr_global_2D[i][j]==8){
          ncp++;
          numcp_global_2D[i][j]=ncp-1;

        }
      }


}
