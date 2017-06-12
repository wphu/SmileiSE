#include "PusherBoris_GC0.h"

#include <iostream>
#include <cmath>

#include "Particles.h"

using namespace std;

PusherBoris_GC0::PusherBoris_GC0(PicParams& params, int ispec)
    : Pusher(params, ispec)
{
}

PusherBoris_GC0::~PusherBoris_GC0()
{
}


/***********************************************************************
	Lorentz Force -- leap-frog (Boris) GC0 scheme -- non-relativistic  ---by wphu
    !!! Now the method is only for the case: 1D3V, the magnetic field line parallels to the x direction
***********************************************************************/
void PusherBoris_GC0::operator() (Particles &particles, int ipart, LocalFields Epart, LocalFields Bpart)
{
    particles.momentum(0, ipart) += charge_over_mass_*Epart.x*dt;

    // Move the particle
    for ( int i = 0 ; i<nDim_ ; i++ ) {
        particles.position_old(i, ipart)  = particles.position(i, ipart);
        particles.position(i, ipart)     += dt*particles.momentum(i, ipart);
    }

}
