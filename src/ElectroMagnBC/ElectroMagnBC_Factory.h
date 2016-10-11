#ifndef ELECTROMAGNBC_FACTORY_H
#define ELECTROMAGNBC_FACTORY_H

#include <vector>

#include "ElectroMagnBC.h"
#include "ElectroMagnBC1D_EC.h"
#include "PicParams.h"

using namespace std;

class ElectroMagnBCFactory {
public:
	static vector<ElectroMagnBC*> create(PicParams& params) {
		vector<ElectroMagnBC*> emBoundCond;

		if ( params.geometry == "1d3v" && params.bcType == "ExternCircuit" ) {
			emBoundCond.resize(1, NULL);
			emBoundCond[0] = new ElectroMagnBC1D_EC(params);
		}
		return emBoundCond;
	}
};


#endif
