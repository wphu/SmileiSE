#include "PusherBoris_imp.h"

#include <iostream>
#include <cmath>

#include "Particles.h"

using namespace std;

PusherBoris_imp::PusherBoris_imp(PicParams& params, int ispec)
    : Pusher(params, ispec)
{
    double Omega_square_p1;
    theta = 1.0;
    factor0 = theta / 2.0;
    factor1 = 1.0 - factor0;

    Omega[0] = 0.5 * charge_over_mass_ * dt * params.externB[0];
    Omega[1] = 0.5 * charge_over_mass_ * dt * params.externB[1];
    Omega[2] = 0.5 * charge_over_mass_ * dt * params.externB[2];
    Omega_square_p1 = 1.0 + pow(Omega[0], 2) + pow(Omega[1], 2) + pow(Omega[2], 2);
    T[0][0] = ( 1.0 + pow(Omega[0], 2) )        / Omega_square_p1;
    T[0][1] = ( Omega[0]*Omega[1] + Omega[2] )  / Omega_square_p1;
    T[0][2] = ( Omega[0]*Omega[2] - Omega[1] )  / Omega_square_p1;
    T[1][0] = ( Omega[0]*Omega[1] - Omega[2] )  / Omega_square_p1;
    T[1][1] = ( 1.0 + pow(Omega[1], 2) )        / Omega_square_p1;
    T[1][2] = ( Omega[1]*Omega[2] + Omega[0] )  / Omega_square_p1;
    T[2][0] = ( Omega[0]*Omega[2] + Omega[1] )  / Omega_square_p1;
    T[2][1] = ( Omega[1]*Omega[2] - Omega[0] )  / Omega_square_p1;
    T[2][2] = ( 1.0 + pow(Omega[2], 2) )        / Omega_square_p1;

}

PusherBoris_imp::~PusherBoris_imp()
{
}



/***********************************************************************
	First push -- leap-frog (Boris) implicit scheme -- non-relativistic  ---by wphu
***********************************************************************/
void PusherBoris_imp::firstPush (Particles &particles, int ipart, LocalFields Epart)
{
    // --------------------------------------
    // SOLVE THE PARTICLE EQUATION OF MOTIONS
    // --------------------------------------

    v_n[0] = particles.momentum(0, ipart);
    v_n[1] = particles.momentum(1, ipart);
    v_n[2] = particles.momentum(2, ipart);

    Au_imp[0] = particles.au_imp(0, ipart);
    Au_imp[1] = particles.au_imp(1, ipart);
    Au_imp[2] = particles.au_imp(2, ipart);

    // no ET term in equation (6)
    S1[0] = v_n[0] + 0.5 * Au_imp[0] * dt + v_n[1] * Omega[2] - v_n[2] * Omega[1];
    S1[1] = v_n[1] + 0.5 * Au_imp[1] * dt + v_n[2] * Omega[0] - v_n[0] * Omega[2];
    S1[2] = v_n[2] + 0.5 * Au_imp[2] * dt + v_n[0] * Omega[1] - v_n[1] * Omega[0];

    v_np1[0] = T[0][0] * S1[0] + T[0][1] * S1[1] + T[0][2] * S1[2];
    v_np1[1] = T[1][0] * S1[0] + T[1][1] * S1[1] + T[1][2] * S1[2];
    v_np1[2] = T[2][0] * S1[0] + T[2][1] * S1[1] + T[2][2] * S1[2];

    particles.momentum(0, ipart) = v_np1[0];
    particles.momentum(1, ipart) = v_np1[1];
    particles.momentum(2, ipart) = v_np1[2];

    // Move the particle
    for ( int i = 0 ; i<nDim_ ; i++ ) {
        particles.position_old(i, ipart)  = particles.position(i, ipart);
        particles.position(i, ipart)     += v_np1[i] * dt;
    }

}



/***********************************************************************
	Second push -- leap-frog (Boris) implicit scheme -- non-relativistic  ---by wphu
***********************************************************************/
void PusherBoris_imp::secondPush (Particles &particles, int ipart, LocalFields Epart)
{
    // --------------------------------------
    // SOLVE THE PARTICLE EQUATION OF MOTIONS
    // --------------------------------------

    v_n[0] = particles.momentum(0, ipart);
    v_n[1] = particles.momentum(1, ipart);
    v_n[2] = particles.momentum(2, ipart);

    Al_imp_np1[0] = charge_over_mass_ * Epart.x;
    Al_imp_np1[1] = charge_over_mass_ * Epart.y;
    Al_imp_np1[2] = charge_over_mass_ * Epart.z;

    S2[0] = Al_imp_np1[0] * dts2;
    S2[1] = Al_imp_np1[1] * dts2;
    S2[2] = Al_imp_np1[2] * dts2;

    deltav[0] = T[0][0] * S2[0] + T[0][1] * S2[1] + T[0][2] * S2[2];
    deltav[1] = T[1][0] * S2[0] + T[1][1] * S2[1] + T[1][2] * S2[2];
    deltav[2] = T[2][0] * S2[0] + T[2][1] * S2[1] + T[2][2] * S2[2];

    v_np1[0] = v_n[0] + deltav[0];
    v_np1[1] = v_n[1] + deltav[1];
    v_np1[2] = v_n[2] + deltav[2];

    particles.momentum(0, ipart) = v_np1[0];
    particles.momentum(1, ipart) = v_np1[1];
    particles.momentum(2, ipart) = v_np1[2];

    // Move the particle
    for ( int i = 0 ; i<nDim_ ; i++ ) {
        particles.position_old(i, ipart)  = particles.position(i, ipart);
        particles.position(i, ipart)     += deltav[i] * dt;

        particles.au_imp(i, ipart) = factor0 * Al_imp_np1[i] + factor1 * particles.al_imp(i, ipart);
        particles.al_imp(i, ipart) = factor1 * Al_imp_np1[i] + factor0 * particles.al_imp(i, ipart);
    }

}
