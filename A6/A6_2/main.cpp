// 04/16/2014 - ENGR 2450 - Meine, Joel
// Lorenz Equations (Pg. 816 | Chapra & Canale), Problem 28.19

// Runge-Kutta Method (Fourth Order) - System of Equations

#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
using namespace std;

const double o = 10, b = 2.666667, r = 28;

vector<double> Lorenz(double t,vector<double> y)
{
	double dxdt,dydt,dzdt;
	dxdt = -o*y[0] + o*y[1];
	dydt = r*y[0] - y[1] - y[0]*y[2];
	dzdt = -b*y[2] + y[0]*y[1];
	vector<double> F = {dxdt,dydt,dzdt};
	return(F);
}

const double f = 60, L = 30, E = 1.25e8, I = 0.05;

vector<double> Sailboat(double z,vector<double> y)
{
	double dy1dz,dy2dz;
	dy1dz = y[1];
	dy2dz = (f/(2*E*I))*pow(L-z,2);
	vector<double> F = {dy1dz,dy2dz};
	return(F);
}

vector<double> Derivs(double x,vector<double> y,int p)
{
	if (p == 1)
	{
		return(Lorenz(x,y));
	}
	if (p == 2)
	{
		return(Sailboat(x,y));
	}
}

void RK4(double& x,vector<double>& y,int& n,double& h,int& p)
{
	vector<double> k1,k2,k3,k4;
	vector<double> ym1,ym2,ye,slope;
	k1 = Derivs(x,y,p);
	for (int i = 0; i < n; i++)
	{
		ym1.push_back(y[i] + k1[i]*(h/2.0));
	}
	k2 = Derivs(x+(h/2.0),ym1,p);
	for (int i = 0; i < n; i++)
	{

		ym2.push_back(y[i] + k2[i]*(h/2.0));
	}
	k3 = Derivs(x+(h/2.0),ym2,p);
	for (int i = 0; i < n; i++)
	{
		ye.push_back(y[i] + k3[i]*h);
	}
	k4 = Derivs(x+h,ye,p);
	for (int i = 0; i < n; i++)
	{
		slope.push_back((k1[i] + 2.0*(k2[i] + k3[i]) + k4[i])/6.0);
		y[i] = y[i] + slope[i]*h;
	}
	x = x + h;
}

void Integrator(double& x,vector<double>& y,int& n,double& h,double& xend,int p)
{
	do {
		if (xend - x < h) h = xend - x;
		RK4(x,y,n,h,p);
	} while (x < xend);
}

int main()
{
	// Lorenz Equations (Pg. 816 | Chapra & Canale)

	vector<double> xp1;
	vector< vector<double> > yp1;
	vector<double> y1,yi1;
	int n1;
	double x1,xi1,xf1,dx1,h1,xout1,xend1;

	yi1 = {5,5,5}; n1 = 3; xi1 = 0; xf1 = 20; dx1 = 0.1; xout1 = 2.0;
	const int size1 = 11;

	x1 = xi1;
	xp1.push_back(x1);
	for (int i = 0; i < n1; i++)
	{
		vector<double> ypw1;
		ypw1.push_back(yi1[i]);
		yp1.push_back(ypw1);
		y1.push_back(yi1[i]);
	}
	do {
		xend1 = x1 + xout1;
		if (xend1 > xf1) xend1 = xf1;
		h1 = dx1;
		Integrator(x1,y1,n1,h1,xend1,1);
		xp1.push_back(x1);
		for (int i = 0; i < n1; i++)
		{
			yp1[i].push_back(y1[i]);
		}
	} while (x1 < xf1);

	std::cout << "Lorenz Equations (Pg. 816 | Chapra & Canale)" << std::endl;
	std::cout << "=================================================" << std::endl;
	std::cout << "dx/dt = -o*x + o*y" << std::endl;
	std::cout << std::endl;
	std::cout << "dy/dt = r*x - y - x*z" << std::endl;
	std::cout << std::endl;
	std::cout << "dz/dt = -b*z + x*y" << std::endl;
	std::cout << std::endl;
	std::cout << "o = " << o << "; b = " << b << "; r = " << r <<std::endl;
	std::cout << std::endl;
	std::cout << "Atmospheric Fluid Motion, x" << std::endl;
	std::cout << "Temperature Variation (Horizontal Axis), y" << std::endl;
	std::cout << "Temperature Variation (Vertical Axis), z" << std::endl;
	std::cout << "Time, t" << std::endl;
	std::cout << "*************************************************" << std::endl;
	std::cout << " t   x         y          z" << std::endl;
	std::cout << "-------------------------------------------------" << std::endl;
	for (int i = 0; i < size1; i++)
	{
		cout << setw(3) << xp1[i];
		cout << setw(9) << setprecision(5) << yp1[0][i];
		cout << setw(11) << setprecision(5) << yp1[1][i];
		cout << setw(9) << setprecision(5) << yp1[2][i] << endl;
	}
	std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	cout << "\n";

	// Problem 28.19

	vector<double> xp2;
	vector< vector<double> > yp2;
	vector<double> y2,yi2;
	int n2;
	double x2,xi2,xf2,dx2,h2,xout2,xend2;

	yi2 = {0,0}; n2 = 2; xi2 = 0; xf2 = L; dx2 = 0.5; xout2 = 5.0;
	const int size2 = 7;

	x2 = xi2;
	xp2.push_back(x2);
	for (int i = 0; i < n2; i++)
	{
		vector<double> ypw2;
		ypw2.push_back(yi2[i]);
		yp2.push_back(ypw2);
		y2.push_back(yi2[i]);
	}
	do {
		xend2 = x2 + xout2;
		if (xend2 > xf2) xend2 = xf2;
		h2 = dx2;
		Integrator(x2,y2,n2,h2,xend2,2);
		xp2.push_back(x2);
		for (int i = 0; i < n2; i++)
		{
			yp2[i].push_back(y2[i]);
		}
	} while (x2 < xf2);

	std::cout << "Problem 28.19" << std::endl;
	std::cout << "=================================================" << std::endl;
	std::cout << "dy1/dz = y2" << std::endl;
	std::cout << std::endl;
	std::cout << "dy2/dz = (f/(2*E*I))*pow(L-z,2)" << std::endl;
	std::cout << std::endl;
	std::cout << "Wind Force, f = " << f << std::endl;
	std::cout << "Modulus of Elasticity, E = " << E << std::endl;
	std::cout << "Mast Length, L = " << L << std::endl;
	std::cout << "Moment of Inertia, I = " << I << std::endl;
	std::cout << std::endl;
	std::cout << "Deflection of Mast (Horizontal Axis), y" << std::endl;
	std::cout << "Deflection of Mast (Vertical Axis), z" << std::endl;
	std::cout << "*************************************************" << std::endl;
	std::cout << " z    y2        dy2/dz" << std::endl;
	std::cout << "-------------------------------------------------" << std::endl;
	for (int i = 0; i < size2; i++)
	{
		cout << setw(3) << xp2[i];
		cout << setw(10) << setprecision(5) << yp2[0][i];
		cout << setw(9) << setprecision(5) << yp2[1][i] << endl;
	}
	std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	cout << "\n";

	system("pause");
	return 0;
}