#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "PSI1D_Sputtering.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
//                                                   MAIN CODE
// ---------------------------------------------------------------------------------------------------------------------
int main (int argc, char* argv[])
{
	double rn[89][1000];
	double z1;
	double z2;
	double am1;
	double am2;
	double es;
	double n;

	double energy, theta, ee, rnion, reion;
	ofstream ofile_rn;

	z1 = 18.0;
	z2 = 74.0;
	am1 = 39.95;
	am2 = 183.85;
	es = 8.7;
	n = 0.11286;

	PSI1D_Sputtering st = PSI1D_Sputtering( z1, am1, z2, am2, es, n );
	for(int i = 1; i < 90; i++)
	{
		for(int j = 1; j < 1001; j++)
		{
			theta = 1.0 * i;
			energy = 1.0 * j;
			rn[i-1][j-1] = st.phy_sput_yield(energy, theta);
		}
	}

	ofile_rn.open("rn.txt");
	for(int i = 1; i < 89; i++)
	{
		for(int j = 1; j < 1001; j++)
		{
			ofile_rn<<setw(20)<<rn[i-1][j-1];
		}
		ofile_rn<<endl;
	}
	ofile_rn.close();
}
