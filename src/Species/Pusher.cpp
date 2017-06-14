#include "Pusher.h"
#include "PicParams.h"

Pusher::Pusher(PicParams& params, int ispec)
{
    mass_          = params.species_param[ispec].mass;
    one_over_mass_ = 1.0/mass_;
    charge_over_mass_ = params.species_param[ispec].charge*one_over_mass_;
    dt             = params.species_param[ispec].timestep_zoom * params.timestep;
    dts2           = params.species_param[ispec].timestep_zoom * params.timestep/2.;

    nDim_          = params.nDim_particle;
    cell_length    = params.cell_length;

}

Pusher::~Pusher()
{
}
