#ifndef ELECTROMAGNBC_H
#define ELECTROMAGNBC_H

#include <string>

using namespace std;

class PicParams;
class SmileiMPI;
class ElectroMagn;
class Solver;

class ElectroMagnBC {
public:
	ElectroMagnBC(){};
	virtual ~ElectroMagnBC(){};

	virtual void apply(ElectroMagn* EMfields, Solver* solver){};

	string bcType;

};


#endif
