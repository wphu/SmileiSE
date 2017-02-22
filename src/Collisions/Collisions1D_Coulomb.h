/*
Collisions1D_Coulomb class
*/

#ifndef COLLISIONS2D_COULOMB_H
#define COLLISIONS2D_COULOMB_H

#include <vector>

#include "Tools.h"
#include "PicParams.h"
#include "InputData.h"
#include "Species.h"
#include "Collisions1D.h"
#include "H5.h"

class Collisions1D_Coulomb : public Collisions1D
{

public:
    //! Constructor for Collisions1D between two species
    Collisions1D_Coulomb(PicParams&,std::vector<Species*>&,SmileiMPI*,unsigned int,std::vector<unsigned int>,std::vector<unsigned int>,double,bool,int);
    ~Collisions1D_Coulomb();


    //! Coulomb logarithm (zero or negative means automatic)
    double coulomb_log;

    //! is true if any of the collisions objects need automatically-computed coulomb log
    bool debye_length_required;
    //! Contains the debye length in each cluster, computed each timestep
    std::vector<double> debye_length_squared;

    //! Method to calculate the Debye length in each cluster
    void calculate_debye_length(PicParams&,std::vector<Species*>&);
    virtual double cos_chi(double);


    virtual void collide_relativistic(PicParams&, SmileiMPI* smpi, std::vector<Species*>&,int);
    // non-relativistic case
    virtual void collide(PicParams&, SmileiMPI* smpi, ElectroMagn* fields, std::vector<Species*>&,int);


private:
    hid_t fileId;

};


#endif
