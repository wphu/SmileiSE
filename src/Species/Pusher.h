#ifndef PUSHER_H
#define PUSHER_H

#include "PicParams.h"
#include "Field.h"

class Particles;


//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class Pusher
{

public:
    //! Creator for Pusher
    Pusher(PicParams& params, int ispec);
    virtual ~Pusher();

    //! Overloading of () operator
    virtual void operator() (Particles &particles, int ipart, LocalFields Epart, LocalFields Bpart, double& gf) = 0;
    virtual void operator() (Particles &particles, int ipart, LocalFields Epart, LocalFields Bpart ) = 0;
    virtual void operator() (Particles &particles, int ipart, LocalFields Epart) = 0;

    // methods for implicit push
    virtual void firstPush (Particles &particles, int ipart, LocalFields Epart) = 0;
    virtual void secondPush (Particles &particles, int ipart, LocalFields Epart) = 0;
protected:
    double dt, dts2;
    //! \todo Move mass_ in Particles_
    // mass_ relative to Species but used in the particle pusher
    double mass_;
    double one_over_mass_;
    double charge_over_mass_;

    int nDim_;

    std::vector<double> cell_length;

};//END class

#endif
