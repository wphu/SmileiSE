/*
// Ref: Wang Hongyu, implicit electrostatic particle in cell/ monte carlo simulation for
//                   the magnetized plasma: Algorithms and application in gas-inductive breakdown
// Equation (6) and (9)

 */

#ifndef PUSHERBORIS_IMP_H
#define PUSHERBORIS_IMP_H

#include "Pusher.h"
#include "userFunctions.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherBoris_imp
//  --------------------------------------------------------------------------------------------------------------------
class PusherBoris_imp : public Pusher {
public:
    //! Creator for Pusher
    PusherBoris_imp(PicParams& params, int ispec);
    ~PusherBoris_imp();

    virtual void firstPush (Particles &particles, int ipart, LocalFields Epart);
    virtual void secondPush (Particles &particles, int ipart, LocalFields Epart);

    // ========== constant ============
    double Omega[3];
    double theta;
    double factor0;
    double factor1;

    // =========== variable ============
    double S1[3];
    double S2[3];
    double deltav[3];
    // acceleration in the n+1 timestep
    double v_n[3];
    double v_np1[3];
    double Al_imp_np1[3];
    double Au_imp[3];


    // =========== tensor =============
    double T[3][3];





};

#endif
