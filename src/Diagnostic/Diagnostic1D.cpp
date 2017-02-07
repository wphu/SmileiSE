#include "Diagnostic1D.h"
#include "Species.h"
#include "SmileiMPI.h"


Diagnostic1D::Diagnostic1D(PicParams& params, SmileiMPI* smpi) :
Diagnostic(params)
{
	particleFlux.resize(n_species);
	heatFlux.resize(n_species);
	angleDist.resize(n_species);
	for(int ispec = 0; ispec < n_species; ispec++)
	{
		particleFlux[ispec].resize(2);
		heatFlux[ispec].resize(2);
		angleDist[ispec].resize(2);
		for(int iDirection = 0; iDirection < 2; iDirection++)
		{
			angleDist[ispec][iDirection].resize(90);
		}
	}
}


void Diagnostic1D::run( SmileiMPI* smpi, vector<Species*>& vecSpecies, int itime )
{
	Species *s1;
	Particles *p1;
	double v_square, v_magnitude;
	int iAngle;
	double flux_temp;
	vector<double> angle_temp;

	angle_temp.resize(90);

	if(itime == dump_step + 1) {
		for(int ispec = 0; ispec < n_species; ispec++)
		{
			for(int iDirection = 0; iDirection < 2; iDirection++)
			{
				particleFlux[ispec][iDirection] = 0.0;
				heatFlux[ispec][iDirection] 	= 0.0;
				for(int iA = 0; iA < 90; iA++)
				{
					angleDist[ispec][iDirection][iA] = 0.0;
				}
			}
		}
	}


	for(int ispec = 0; ispec < n_species; ispec++)
	{
		s1 = vecSpecies[ispec];
		p1 = &(s1->psi_particles);
		for(int iPart = 0; iPart < p1->size(); iPart++)
		{
			v_square = p1->momentum(0,iPart) * p1->momentum(0,iPart) + p1->momentum(1,iPart) * p1->momentum(1,iPart) + p1->momentum(2,iPart) * p1->momentum(2,iPart);
			v_magnitude = sqrt(v_square);
			iAngle = acos( abs(p1->momentum(0,iPart)) / v_magnitude ) / PI_ov_2;
			if( p1->position(0,iPart) < 0.0 ) {
				particleFlux[ispec][0]++;
				heatFlux[ispec][0] += 0.5 * s1->species_param.mass * v_square;
				angleDist[ispec][0][iAngle]++;
			}
			else if( p1->position(0,iPart) > sim_length[0] ) {
				particleFlux[ispec][1]++;
				heatFlux[ispec][1] += 0.5 * s1->species_param.mass * v_square;
				angleDist[ispec][1][iAngle]++;
			}
		}

		if(itime == dump_step) {
			particleFlux[ispec][0] /= (timestep * dump_step);
			particleFlux[ispec][1] /= (timestep * dump_step);
			heatFlux[ispec][0] /= (timestep * dump_step);
			heatFlux[ispec][1] /= (timestep * dump_step);
		}

	}

	if(itime == dump_step) {
		for(int ispec = 0; ispec < n_species; ispec++)
		{
			for(int iDirection = 0; iDirection < 2; iDirection++)
			{
				MPI_Allreduce( &particleFlux[ispec][iDirection], &flux_temp, 1, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
				particleFlux[ispec][iDirection] = flux_temp;

				MPI_Allreduce( &heatFlux[ispec][iDirection], &flux_temp, 1, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
				heatFlux[ispec][iDirection] = flux_temp;

				MPI_Allreduce( &angleDist[ispec][iDirection][0], &angle_temp[0], 90, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
				angleDist[ispec][iDirection] = angle_temp;
			}
		}
	}



}
