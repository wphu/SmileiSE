#ifndef ARRAY4D_H
#define ARRAY4D_H

#include <cmath>

#include <vector>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>

#include "Tools.h"

//! Class Array4D: generic class allowing to define vectors
class Array4D
{

public:

    //! name of the Array4D
    std::string name;

    //! Constructor for Array4D: with no input argument
    Array4D() {
        ;
    };

    //! Constructor for Array4D: with the Array4D dimensions as input argument
    Array4D( std::vector<unsigned int> dims ) {
        ;
    };
    //! Constructor, isPrimal define if mainDim is Primal or Dual
    Array4D( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal ) {
        ;
    };

    //! Constructor for Array4D: with the Array4D dimensions and dump file name as input argument
    Array4D( std::vector<unsigned int> dims, std::string name_in ) : name(name_in) {
        ;
    } ;

    //! Constructor for Array4D: isPrimal define if mainDim is Primal or Dual
    Array4D( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal, std::string name_in ) : name(name_in) {
        ;
    } ;

    //! Destructor for Array4D
    ~Array4D() {
        ;
    } ;

    //! Virtual method used to allocate Array4D
    void allocateDims(std::vector<unsigned int> dims)
    {
        dims_=dims;
        if (dims_.size()!=4) ERROR("Alloc error must be 4 : " << dims.size());
        if (data_) delete [] data_;

        data_ = new double[dims_[0]*dims_[1]*dims_[2]*dims_[3]];
        //! \todo{check row major order!!!}
        data_4D= new double***[dims_[0]*dims_[1]*dims_[2]];
        for (unsigned int i=0; i<dims_[0]; i++)
        {
            data_4D[i]= new double**[dims_[1]];
            for (unsigned int j=0; j<dims_[1]; j++)
            {
                data_4D[i][j] = new double*[dims_[2]];
                for (unsigned int k=0; k<dims_[2]; k++)
                {
                    data_4D[i][j][k] = data_ + i*dims_[1]*dims_[2]*dims_[3] + j*dims_[2]*dims_[3] + k*dims_[3];
                }
            }
        }//i

        //DEBUG(10,"Fields 3D created: " << dims_[0] << "x" << dims_[1] << "x" << dims_[2]);
        globalDims_ = dims_[0]*dims_[1]*dims_[2]*dims_[3];
    };


    //! vector containing the dimensions of the Array4D
    //! \todo private/friend/modify SmileiMPI* (JD)
    std::vector<unsigned int> dims_;


    //! returns the dimension of the Array4D
	inline std::vector<unsigned int> dims () {return dims_;}
    //! All arrays may be viewed as a 1D array
    //! Linearized diags
    unsigned int globalDims_;
    //! pointer to the linearized array
    double* data_;

    // data_4D[nx][ny][nz][nv]
    double ****data_4D;

    //! method used to put all entry of a Array4D at a given value val
    inline void put_to(double val)
    {
        if (data_)
            for (unsigned int i=0; i<globalDims_; i++) data_[i] = val;
    }


    //! 2D reference access to the linearized array (with check in DEBUG mode)
    inline double& operator () (unsigned int i,unsigned int j,unsigned int k,unsigned int l)
    {
        return data_4D[i][j][k][l];
    };
    //! 2D access to the linearized array (with check in DEBUG mode)
    inline double operator () (unsigned int i, unsigned int j,unsigned int k,unsigned int l) const
    {
        return data_4D[i][j][k][l];
    };


protected:

private:

};

#endif
