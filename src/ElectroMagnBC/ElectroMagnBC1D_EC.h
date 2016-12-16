#ifndef ELECTROMAGNBC1D_EC_H
#define ELECTROMAGNBC1D_EC_H

#include "ElectroMagnBC.h"

class PicParams;
class ElectroMagn;
class Solver;

// EC: Extern Circuit
class ElectroMagnBC1D_EC : public ElectroMagnBC {
public:
	ElectroMagnBC1D_EC(PicParams&);
	~ElectroMagnBC1D_EC(){};

	virtual void apply(ElectroMagn* EMfields, Solver* solver);

	// capacitance between the discharge gap
	double capacitance;
};


#endif
