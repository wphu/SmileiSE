#include "ElectroMagnBC1D_EC.h"

#include "PicParams.h"
#include "ElectroMagn.h"
#include "SolverFactory.h"

ElectroMagnBC1D_EC::ElectroMagnBC1D_EC(PicParams& params)
{
	bcType = params.bcType;
}


void ElectroMagnBC1D_EC::apply(ElectroMagn* EMfields, Solver* solver)
{
	double voltage_right;
	if ( bcType == "ExternCircuit" ) {
		voltage_right = ( EMfields->totCharge[1][0] + EMfields->depCharge[1][0] + EMfields->emitCharge[1][0]) /
						capacitance;
	}
	solver->bc_e_value[0][1] = voltage_right;
}
