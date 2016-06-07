/*
Collisions2D class
*/

#ifndef COLLISIONS2D_H
#define COLLISIONS2D_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "H5.h"
#include "Collisions.h"


class Collisions2D : public Collisions
{

public:
    //! Constructor for Collisions between two species
    Collisions2D(){};
    virtual ~Collisions2D(){};

protected:
    //! Inverse of the spatial-step
    double dx_inv_;
    double dy_inv_;
    int i_domain_begin;
    int j_domain_begin;

private:
};


#endif
