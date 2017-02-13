#include "Diagnostic.h"

Diagnostic::Diagnostic(PicParams &params) :
n_species(params.species_param.size()),
sim_length(params.sim_length),
dump_step(params.dump_step),
timestep(params.timestep)
{
	PI_ov_2 = 0.5 * params.const_pi;
	const_e = params.const_e;
}
