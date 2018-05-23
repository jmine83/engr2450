// 04/16/2014 - ENGR 2450 - Meine, Joel
// Problem 25.21

// Runge-Kutta Method (Fourth Order)

#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
using namespace std;

vector<double> T = {1950,1960,1970,1980,1990,2000}; // Year
vector<double> P = {2555,3040,3708,4454,5276,6079}; // Population (people in millions)
const int size = 6;

const double kgm = 0.026; // Maximum Growth Rate under Unlimited Conditions
const double pmax = 12000; // Carrying Capacity (people in millions)

double Population(double t,double p)
{
	double dpdt;
	p = pmax / (1-(1-(pmax/P[0]))*exp(-kgm*(t-T[0]))); // Population at Time
	dpdt = kgm*(1-(p/pmax))*p; // Growth Rate of Population with Time
	return(dpdt);
}

double Derivs(double x,double y)
{
	return(Population(x,y));
}

void RK4(double& x,double& y,double& h,double& ynew)
{
	double k1,k2,k3,k4;
	double ym,ye,slope;
	k1 = Derivs(x,y);
	ym = y + k1*(h/2);
	k2 = Derivs(x+(h/2),ym);
	ym = y + k2*(h/2);
	k3 = Derivs(x+(h/2),ym);
	ye = y + k3*h;
	k4 = Derivs(x+h,ye);
	slope = (k1 + 2*(k2 + k3) + k4)/6;
	ynew = y + slope*h;
	x = x + h;
}

void Integrator(double& x,double& y,double& h,double& xend)
{
	double ynew;
	do {
		if (xend - x < h) h = xend - x;
		RK4(x,y,h,ynew);
		y = ynew;
	} while (x < xend);
}

int main()
{
	// Problem 25.21 - Runge-Kutta Method (Fourth Order)

	vector<double> xp,yp;
	double h;
	double x,xi,xf,xend,xout,dx,y;

	xi = T[0]; y = P[0]; xout = 10.0; xf = T.back(); dx = 0.1;

	x = xi;
	xp.push_back(x);
	yp.push_back(y);
	do {
		xend = x + xout;
		if (xend > xf) xend = xf;
		h = dx;
		Integrator(x,y,h,xend);
		xp.push_back(x);
		yp.push_back(y);
	} while (x < xf);

	std::cout << "Chapter 25 - Problem 25.21" << std::endl;
	std::cout << "===================================================================" << std::endl;
	std::cout << "Growth Rate of Population with Time (millions of people per year)" << std::endl;
	std::cout << std::endl;
	std::cout << "dp/dt = kgm * (1-(p/pmax))*p" << std::endl;
	std::cout << std::endl;
	std::cout << "Population at Time (millions of people)" << std::endl;
	std::cout << std::endl;
	std::cout << "p = pmax / (1-(1-(pmax/p0))*exp(-kgm*(t-t0)))" << std::endl;
	std::cout << std::endl;
	std::cout << "Actual Population at Time (people in millions), pa" << std::endl;
	std::cout << "Initial Population at t0 (people in millions), p0 = " << P[0] << std::endl;
	std::cout << "Initial Year with p0, t0 = " << T[0] << std::endl;
	std::cout << "Maximum Growth Rate under Unlimited Conditions, kgm = " << kgm << std::endl;
	std::cout << "Carrying Capacity (people in millions), pmax = " << pmax << std::endl;
	std::cout << "*******************************************************************" << std::endl;
	std::cout << " t     pa    p" << std::endl;
	std::cout << "-------------------------------------------------------------------" << std::endl;
	for (int i = 0; i < size; i++)
	{
		cout << setw(5) << xp[i];
		cout << setw(6) << P[i];
		cout << setw(11) << setprecision(8) << yp[i] << endl;
	}
	std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	cout << "\n";

	system("pause");
	return 0;
}