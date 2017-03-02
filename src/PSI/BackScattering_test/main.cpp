#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "PSI1D_Backscattering.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
//                                                   MAIN CODE
// ---------------------------------------------------------------------------------------------------------------------
int main (int argc, char* argv[])
{
	double rn[91][500], re[91][500];
	int nz1;         // atomic number of incident atomic
	int m1;          // atomic mass of incident atomic (amu)
	int ne;             // number of constituent elements in the target.
	vector<int> nz2;	// array for atomic numbers of the constituents.
	vector<int> nw;     // array for relative numbers of the constituents.
	double energy, theta, ee, rnion, reion;
	ofstream ofile_rn, ofile_re;

	nz1 = 20;
	m1 = 40;
	ne = 1;
	nz2.push_back(74);
	nw.push_back(1);

	PSI1D_Backscattering bs = PSI1D_Backscattering( nz1, m1, ne, nz2, nw );
	for(int i = 0; i < 91; i++)
	{
		for(int j = 0; j < 500; j++)
		{
			theta = 1.0 * i;
			energy = 0.1 * (j+1);
			bs.scatter(rnion, reion, theta, energy);
			rn[i][j] = rnion;
			re[i][j] = reion;
		}
	}

	ofile_rn.open("rn.txt");
	ofile_re.open("re.txt");
	for(int i = 0; i < 91; i++)
	{
		for(int j = 0; j < 500; j++)
		{
			ofile_rn<<setw(10)<<rn[i][j];
			ofile_re<<setw(10)<<re[i][j];
		}
		ofile_rn<<endl;
		ofile_re<<endl;
	}
	ofile_rn.close();
	ofile_re.close();
}
