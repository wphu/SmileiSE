/*! @file Pusher.h

 @brief Pusher.h  generic class for the particle pusher

 @date 2013-02-15
 */

#ifndef PUSHERBORIS_H
#define PUSHERBORIS_H

#include "Pusher.h"
#include "userFunctions.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherBoris
//  --------------------------------------------------------------------------------------------------------------------
class PusherBoris : public Pusher {
public:
    //! Creator for Pusher
    PusherBoris(PicParams& params, int ispec);
    ~PusherBoris();
    //! Overloading of () operator
    virtual void operator() (Particles &particles, int ipart, LocalFields Epart, LocalFields Bpart, double& gf);
    virtual void operator() (Particles &particles, int ipart, LocalFields Epart, LocalFields Bpart);
    virtual void operator() (Particles &particles, int ipart, LocalFields Epart);

    // methods for implicit push
    virtual void firstPush (Particles &particles, int ipart, LocalFields Epart){};
    virtual void secondPush (Particles &particles, int ipart, LocalFields Epart){};
};

#endif
