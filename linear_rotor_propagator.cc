#include <iostream>
#include <stdlib.h>
#include<cmath>
#include <fstream>

using namespace std;

extern "C" 
{
	double plgndr(int j,int m,double x);
}

int main(int argc, char *argv[])
{
	double temperature, bconst;
	int nslice, size_theta, size_theta1, size_phi1;
	string spin_isomer;

	if (argc != 8)
	{
		cout<<"Command line arguments mismatch!"<<endl;
		exit(1);
	}
	else
	{	
		temperature = atof(argv[1]);
		nslice      = atoi(argv[2]);
		bconst      = atof(argv[3]);
		size_theta  = atoi(argv[4]);
		spin_isomer = argv[5];
		size_theta1 = atoi(argv[6]);
		size_phi1   = atoi(argv[7]);
	}

	string isomer, basis_type;
	if (spin_isomer == "spinless") {
		isomer = "-";
		basis_type = "";
	}
	if (spin_isomer == "para") {
		isomer = "-p-";
		basis_type = "even";
	}
	if (spin_isomer == "ortho") {
		isomer = "-o-";
		basis_type = "odd";
	}

	double CMRECIP2KL = 1.4387672;       	//cm^-1 to Kelvin conversion factor
	double beta = 1.0/temperature;
	double tau = beta/nslice;
	bconst = bconst*CMRECIP2KL;

	double tol = 10e-16;
	int maxj=1000;

	int id=0;
	int jmax;

	double exp1;
	for (int j=0; j<maxj; j++) 
	{ 
		exp1 = exp(-tau*bconst*j*(j+1.0));
		if (exp1 < tol) 
		{
			cout<<j<<", "<<exp1<<", "<<-tau*bconst*j*(j+1.0)<<endl;
			jmax=j;
			cout<<"jmax = "<<jmax<<endl;
			break;
		}
	}
	if (id == (maxj-1)) 
	{
		jmax=maxj;
		cout<<"!!! Warning: maxj is reached"<<endl;
	}

	double cost[size_theta];
	double dens[size_theta];
	double erot[size_theta];
	double erotsq[size_theta];

	ofstream fid;
	fid.open("linden_cc.out");
	if( !fid ) { 
		cerr << "Error: file could not be opened" << endl;
		exit(1);
	}

	string str_file = "#First col. --> ei.ej; 2nd, 3rd and 4th columns are the density and energy and heat capacity estimators, respectively.";
	fid<<str_file<<endl;

	double cstep=2.0/(double)(size_theta-1);
	double Pn, tmp, sum, sum_eng, sum_engsq;
	int Nj;
	for (int ic=0; ic<size_theta; ic++) 
	{
		cost[ic] = ic*cstep-1.0;
		sum = 0.0;
		sum_eng = 0.0;
		sum_engsq = 0.0;
		for (int j=0; j<(jmax+1); j++)
		{
			Nj = (2.0*j+1.0);
			Pn = plgndr(j,0,cost[ic]);
			tmp = Nj*Pn*exp(-tau*bconst*j*(j+1.0));
			sum += tmp;
			sum_eng += tmp*j*(j+1.0);
			sum_engsq +=tmp*j*(j+1.0)*j*(j+1.0);
		}
		dens[ic] = sum/(4.0*M_PI);
		erot[ic] = sum_eng/(4.0*M_PI*dens[ic]*nslice);
		erotsq[ic] = sum_engsq/(4.0*M_PI*nslice*nslice*dens[ic]);
		fid<<cost[ic]<<"  "<<dens[ic]<<"  "<<erot[ic]<<"  "<<erotsq[ic]<<endl;
	}
	fid.close();


	exit(1);
}


