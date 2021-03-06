#ifndef DIAGNOSTIC1D_H
#define DIAGNOSTIC1D_H

#include "Diagnostic.h"
#include "PicParams.h"
#include "SmileiMPI.h"
#include "Field1D.h"

class Diagnostic1D : public Diagnostic {

public :

    Diagnostic1D(PicParams& params, SmileiMPI* smpi, ElectroMagn* EMfields);
    virtual ~Diagnostic1D();

    //! Runs the diag for all patches for local diags.
    virtual void run( SmileiMPI* smpi, vector<Species*>& vecSpecies, ElectroMagn* EMfields, int itime ) ;
    // calculate velocity and temperature of each species
	 void calVT(SmileiMPI* smpi, vector<Species*>& vecSpecies, ElectroMagn* EMfields, int itime);

   // calculate total energy(particles and electric field)
   void calTotalEnergy(SmileiMPI* smpi, vector<Species*>& vecSpecies, ElectroMagn* EMfields, int itime);


	vector< vector<double> > particleFlux;             //particleFlux[ispec][iDirection]
	vector< vector<double> > heatFlux;                 //heatFlux[ispec][iDirection]
    vector< vector< vector<double> > > angleDist;      //angleDist[ispec][iDirection][iAngle]
    vector< int > particleNumber;                      // particleNumber[ispec]
    vector< double > kineticEnergy;                    // kineticEnergy[ispec]
    vector<double> totalParticleEnergy;
    double totalElectricFieldEnergy;



    Field1D *ptclNum1D;


protected :

    //! Inverse of the spatial step 1/dx
    // parameters to project macroscopic velocity and temperature
    double dx_inv_;
    int index_domain_begin;

};

#endif
