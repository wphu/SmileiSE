#include "Diagnostic.h"

Diagnostic::Diagnostic(PicParams &params) :
n_species(params.species_param.size()),
sim_length(params.sim_length),
dump_step(params.dump_step),
timestep(params.timestep)
{
	PI_ov_2 = 0.5 * params.const_pi;
	const_e = params.const_e;

	double B_magnitude = pow(params.externB[0], 2) + pow(params.externB[1], 2) + pow(params.externB[2], 2);
	B_magnitude = sqrt(B_magnitude);
	sinPhi = params.externB[0] / B_magnitude;
	cosPhi = sqrt(1.0 - sinPhi * sinPhi);
}
