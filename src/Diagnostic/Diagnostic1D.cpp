#include "Diagnostic1D.h"
#include "Species.h"
#include "SmileiMPI_Cart1D.h"
#include "ElectroMagn.h"
#include <iomanip>
#include <fstream>

using namespace std;

Diagnostic1D::Diagnostic1D(PicParams& params, SmileiMPI* smpi, ElectroMagn* EMfields) :
Diagnostic(params)
{
	SmileiMPI_Cart1D* smpi1D = static_cast<SmileiMPI_Cart1D*>(smpi);
	dx_inv_  = 1.0/params.cell_length[0];
	index_domain_begin = smpi1D->getCellStartingGlobalIndex(0);

	particleFlux.resize(n_species);
	heatFlux.resize(n_species);
	angleDist.resize(n_species);
	particleNumber.resize(n_species);
	kineticEnergy.resize(n_species);
	totalParticleEnergy.resize(n_species);
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
	ptclNum1D = new Field1D(EMfields->dimPrim, "ptclNum");
}

Diagnostic1D::~Diagnostic1D()
{
	delete ptclNum1D;
}



void Diagnostic1D::run( SmileiMPI* smpi, vector<Species*>& vecSpecies, ElectroMagn* EMfields, int itime )
{
	Species *s1;
	Particles *p1;
	double v_square, v_magnitude;
	double mass_ov_2;
	double wlt;			// weight * cell_length / time for calculating flux
	int iAngle;
	double flux_temp;
	vector<double> angle_temp;

	angle_temp.resize(90);

	// reset diagnostic parameters to zero
	if( ((itime - 1) % dump_step) == 0 ) {
		for(int ispec = 0; ispec < n_species; ispec++)
		{
			for(int iDirection = 0; iDirection < 2; iDirection++)
			{
				particleFlux[ispec][iDirection] = 0.0;
				heatFlux	[ispec][iDirection] = 0.0;
				for(int iA = 0; iA < 90; iA++)
				{
					angleDist[ispec][iDirection][iA] = 0.0;
				}
			}
		}
	}

	// calculate diagnostic parameters
	if( (itime % dump_step) > (dump_step - avg_step) || (itime % dump_step) == 0 )
	{
		for(int ispec = 0; ispec < n_species; ispec++)
		{
			s1 = vecSpecies[ispec];
			p1 = &(s1->psi_particles);

			// 0.5 * mass
			mass_ov_2 = 0.5 * s1->species_param.mass;
			for(int iPart = 0; iPart < p1->size(); iPart++)
			{
				//cout<<"particle number:  "<<p1->position(0,iPart)<<endl;
				v_square = p1->momentum(0,iPart) * p1->momentum(0,iPart) + p1->momentum(1,iPart) * p1->momentum(1,iPart) + p1->momentum(2,iPart) * p1->momentum(2,iPart);
				v_magnitude = sqrt(v_square);
				iAngle = 90.0 * acos( abs(p1->momentum(0,iPart)) / v_magnitude ) / PI_ov_2;
				if( p1->position(0,iPart) < 0.0 ) {
					//cout<<"particle number:  "<<p1->position(0,iPart)<<endl;
					particleFlux[ispec][0]++;
					heatFlux[ispec][0] += mass_ov_2 * v_square;
					if (iAngle >= 0 && iAngle < 90)
					{
						angleDist[ispec][0][iAngle]++;
					}
					else
					{
						WARNING("iAngle out of range: iAngle = "<<iAngle);
					}

				}
				else if( p1->position(0,iPart) > sim_length[0] ) {
					particleFlux[ispec][1]++;
					heatFlux[ispec][1] += mass_ov_2 * v_square;
					if (iAngle >= 0 && iAngle < 90)
					{
						angleDist[ispec][1][iAngle]++;
					}
					else
					{
						WARNING("iAngle out of range: iAngle = "<<iAngle);
					}
				}
			}
		}
	}


	// MPI gather diagnostic parameters to master
	if( (itime % dump_step) == 0 ) {
		for(int ispec = 0; ispec < n_species; ispec++)
		{
			s1 = vecSpecies[ispec];
			wlt = s1->species_param.weight / (dx_inv_ * timestep * avg_step);
			particleFlux[ispec][0] *= wlt;
			particleFlux[ispec][1] *= wlt;
			heatFlux[ispec][0] 	   *= wlt;
			heatFlux[ispec][1]     *= wlt;

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

	// calculate velocity and temperature of each species
	calVT(smpi, vecSpecies, EMfields, itime);

	//calTotalEnergy(smpi, vecSpecies, EMfields, itime);



}


void Diagnostic1D::calVT(SmileiMPI* smpi, vector<Species*>& vecSpecies, ElectroMagn* EMfields, int itime)
{
	Species *s1;
	Particles *p1;
	int i_temp;
	double xjn,xjmxi;
	double m_ov_3e;
	double m_ov_2;
	double vx, vy, vz;
	double avg_step_inv_;
	double v_square;
	vector<int> number_temp;
	vector<double> energy_temp;

	number_temp.resize(n_species);
	energy_temp.resize(n_species);
	avg_step_inv_ = 1.0 / avg_step;
	for(int iSpec = 0; iSpec < vecSpecies.size(); iSpec++)
	{
		s1 = vecSpecies[iSpec];
		p1 = &(s1->particles);
		Field1D* Vx1D_s = static_cast<Field1D*>(EMfields->Vx_s[iSpec]);
		Field1D* Vy1D_s = static_cast<Field1D*>(EMfields->Vy_s[iSpec]);
		Field1D* Vz1D_s = static_cast<Field1D*>(EMfields->Vz_s[iSpec]);
		Field1D* Vp1D_s = static_cast<Field1D*>(EMfields->Vp_s[iSpec]);

		Field1D* Vx1D_s_avg = static_cast<Field1D*>(EMfields->Vx_s_avg[iSpec]);
		Field1D* Vy1D_s_avg = static_cast<Field1D*>(EMfields->Vy_s_avg[iSpec]);
		Field1D* Vz1D_s_avg = static_cast<Field1D*>(EMfields->Vz_s_avg[iSpec]);
		Field1D* Vp1D_s_avg = static_cast<Field1D*>(EMfields->Vp_s_avg[iSpec]);

		Field1D* T1D_s = static_cast<Field1D*>(EMfields->T_s[iSpec]);
		Field1D* T1D_s_avg = static_cast<Field1D*>(EMfields->T_s_avg[iSpec]);

		if( ((itime - 1) % dump_step) == 0 ) {
			Vx1D_s_avg->put_to(0.0);
			Vy1D_s_avg->put_to(0.0);
			Vz1D_s_avg->put_to(0.0);
			Vp1D_s_avg->put_to(0.0);
			T1D_s_avg ->put_to(0.0);
		}

		// calculate macroscopic velocity (average velocity) and particle number at grid points
		if( (itime % dump_step) > (dump_step - avg_step) || (itime % dump_step) == 0 )
		{
			m_ov_3e = s1->species_param.mass / ( const_e * 3.0 );
			m_ov_2 = s1->species_param.mass / 2.0;
			ptclNum1D->put_to(0.0);
			Vx1D_s->put_to(0.0);
			Vy1D_s->put_to(0.0);
			Vz1D_s->put_to(0.0);
			Vp1D_s->put_to(0.0);
			T1D_s->put_to(0.0);

			// reset particleNumber and kineticEnergy
			particleNumber[iSpec] = 0;
			kineticEnergy[iSpec] = 0.0;
			// get particleNumber
			particleNumber[iSpec] = p1->size();

			for(int iPart = 0; iPart < p1->size(); iPart++)
			{
				//Locate particle on the grid
				xjn    = p1->position(0, iPart) * dx_inv_;  	// normalized distance to the first node
				i_temp      = floor(xjn);                   		// index of the central node
				xjmxi  = xjn - (double)i_temp;              		// normalized distance to the nearest grid point

				i_temp -= index_domain_begin;
				(*ptclNum1D)(i_temp) 	+= 1.0;
				(*Vx1D_s)(i_temp) 		+= p1->momentum(0, iPart);
				(*Vy1D_s)(i_temp) 		+= p1->momentum(1, iPart);
				(*Vz1D_s)(i_temp) 		+= p1->momentum(2, iPart);
				(*Vp1D_s)(i_temp) 		+= (p1->momentum(0, iPart) * sinPhi + p1->momentum(1, iPart) * cosPhi);

				// calculate total kineticEnergy
				v_square = p1->momentum(0, iPart) * p1->momentum(0, iPart) + p1->momentum(1, iPart) * p1->momentum(1, iPart) + p1->momentum(2, iPart) * p1->momentum(2, iPart);
				kineticEnergy[iSpec] += ( m_ov_2 * v_square );
			}
			for(int i = 0; i < ptclNum1D->dims_[0]; i++)
			{
				if( (*ptclNum1D)(i) != 0.0 )
				{
					(*Vx1D_s)(i) /= (*ptclNum1D)(i);
					(*Vy1D_s)(i) /= (*ptclNum1D)(i);
					(*Vz1D_s)(i) /= (*ptclNum1D)(i);
					(*Vp1D_s)(i) /= (*ptclNum1D)(i);
				}
			}

			// calculate temperature
			for(int iPart = 0; iPart < p1->size(); iPart++)
			{
				//Locate particle on the grid
				xjn    = p1->position(0, iPart) * dx_inv_;  	// normalized distance to the first node
				i_temp      = floor(xjn);                   		// index of the central node
				xjmxi  = xjn - (double)i_temp;              		// normalized distance to the nearest grid point

				i_temp -= index_domain_begin;
				//vx = p1->momentum(0, iPart);
				//vy = p1->momentum(1, iPart);
				//vz = p1->momentum(2, iPart);
				vx = p1->momentum(0, iPart) - (*Vx1D_s)(i_temp);
				vy = p1->momentum(1, iPart) - (*Vy1D_s)(i_temp);
				vz = p1->momentum(2, iPart) - (*Vz1D_s)(i_temp);
				(*T1D_s)(i_temp) 		+= ( vx * vx + vy * vy + vz * vz );
			}
			for(int i = 0; i < ptclNum1D->dims_[0]; i++)
			{
				if( (*ptclNum1D)(i) != 0.0 )
				{
					(*T1D_s)(i) = (*T1D_s)(i) * m_ov_3e / (*ptclNum1D)(i);
				}

			}

			// sum velocity and temperature
			for(int i = 0; i < ptclNum1D->dims_[0]; i++)
			{
				(*Vx1D_s_avg)(i) += (*Vx1D_s)(i);
				(*Vy1D_s_avg)(i) += (*Vy1D_s)(i);
				(*Vz1D_s_avg)(i) += (*Vz1D_s)(i);
				(*Vp1D_s_avg)(i) += (*Vp1D_s)(i);
				(*T1D_s_avg)(i)  += (*T1D_s)(i);
			}
		}


		// Calculate the average parameters and MPI gather
		if( (itime % dump_step) == 0 )
		{
			for(int i = 0; i < ptclNum1D->dims_[0]; i++)
			{
				(*Vx1D_s_avg)(i) *= avg_step_inv_;
				(*Vy1D_s_avg)(i) *= avg_step_inv_;
				(*Vz1D_s_avg)(i) *= avg_step_inv_;
				(*Vp1D_s_avg)(i) *= avg_step_inv_;
				(*T1D_s_avg)(i)  *= avg_step_inv_;
			}

			// another way: firstly gather V, T, ptclNum, then calculate V_global, T_global
			SmileiMPI_Cart1D* smpi1D = static_cast<SmileiMPI_Cart1D*>(smpi);
			smpi1D->gatherRho( static_cast<Field1D*>(EMfields->Vx_s_global_avg[iSpec]), Vx1D_s_avg );
			smpi1D->gatherRho( static_cast<Field1D*>(EMfields->Vy_s_global_avg[iSpec]), Vy1D_s_avg );
			smpi1D->gatherRho( static_cast<Field1D*>(EMfields->Vz_s_global_avg[iSpec]), Vz1D_s_avg );
			smpi1D->gatherRho( static_cast<Field1D*>(EMfields->Vp_s_global_avg[iSpec]), Vp1D_s_avg );
			smpi1D->gatherRho( static_cast<Field1D*>(EMfields->T_s_global_avg [iSpec]), T1D_s_avg );

		}

	}
	// gather particleNumber and kineticEnergy to master process
	if( (itime % dump_step) == 0 )
	{
		MPI_Allreduce( &particleNumber[0], &number_temp[0], n_species , MPI_INT,MPI_SUM, MPI_COMM_WORLD);
		particleNumber = number_temp;
		MPI_Allreduce( &kineticEnergy[0], &energy_temp[0], n_species , MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
		kineticEnergy = energy_temp;
	}


}


void Diagnostic1D::calTotalEnergy(SmileiMPI* smpi, vector<Species*>& vecSpecies, ElectroMagn* EMfields, int itime)
{
		Species *s1;
		Particles *p1;
		double m_ov_2, v_square;
		int i;
		double totalElectricFieldEnergy_temp;
		vector<double> totalParticleEnergy_temp;
		totalParticleEnergy_temp.resize(n_species);

		for(int iSpec = 0; iSpec < vecSpecies.size(); iSpec++)
		{
				s1 = vecSpecies[iSpec];
				p1 = &(s1->particles);

				m_ov_2 = s1->species_param.mass / 2.0;

				// reset particleNumber and kineticEnergy
				totalParticleEnergy[iSpec] = 0.0;

				for(int iPart = 0; iPart < p1->size(); iPart++)
				{
					// calculate total kineticEnergy
					v_square = p1->momentum(0, iPart) * p1->momentum(0, iPart) + p1->momentum(1, iPart) * p1->momentum(1, iPart) + p1->momentum(2, iPart) * p1->momentum(2, iPart);
					totalParticleEnergy[iSpec] += ( m_ov_2 * v_square );
				}

		}

		totalElectricFieldEnergy = 0.0;
		Field1D* Ex1D = static_cast<Field1D*>(EMfields->Ex_);
		for(int i = oversize[0]; i < Ex1D->dims_[0] - oversize[0]; i++)
		{
				totalElectricFieldEnergy += (*Ex1D)(i) * (*Ex1D)(i);
		}

		smpi->reduceDoubleVector(&totalParticleEnergy[0], &totalParticleEnergy_temp[0], n_species);
		smpi->reduceDoubleVector(&totalElectricFieldEnergy, &totalElectricFieldEnergy_temp, 1);
		totalParticleEnergy = totalParticleEnergy_temp;
		totalElectricFieldEnergy = totalElectricFieldEnergy_temp;

		if(smpi->isMaster())
		{
				ofstream outfile;
				outfile.open("totalEnergy.txt", ios::app);
				outfile<<setw(15)<<itime<<setw(15)<<totalElectricFieldEnergy<<setw(15)<<totalParticleEnergy[0]<<setw(15)<<totalParticleEnergy[1]
							 <<setw(15)<<totalParticleEnergy[0]+totalParticleEnergy[1]<<endl;
				outfile.close();
		}

}
