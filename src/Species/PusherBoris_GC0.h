/*
Guiding Center pusher 0: the charge is all located in the guiding center
This is usually appropriate for the case: radius of gyration is much less than the grid cell length
*/

#ifndef PUSHERBORIS_GC0_H
#define PUSHERBORIS_GC0_H

#include "Pusher.h"
#include "userFunctions.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherBoris
//  --------------------------------------------------------------------------------------------------------------------
class PusherBoris_GC0 : public Pusher {
public:
    //! Creator for Pusher
    PusherBoris_GC0(PicParams& params, int ispec);
    ~PusherBoris_GC0();
    //! Overloading of () operator
    virtual void operator() (Particles &particles, int ipart, LocalFields Epart, LocalFields Bpart, double& gf){};
    virtual void operator() (Particles &particles, int ipart, LocalFields Epart, LocalFields Bpart);
    virtual void operator() (Particles &particles, int ipart, LocalFields Epart){};

    // methods for implicit push
    virtual void firstPush (Particles &particles, int ipart, LocalFields Epart){};
    virtual void secondPush (Particles &particles, int ipart, LocalFields Epart){};
};

#endif
