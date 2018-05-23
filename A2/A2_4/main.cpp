// 02/05/2014 - ENGR 2450 - Meine, Joel
// Chapter 6 - Problem 6.18

#include <iostream>
#include <math.h>
using namespace std;

const double V = 1 * pow(10,6); // Volume (m^3)
const double Q = 1 * pow(10,5); // Volume per Year (m^3/yr)
const double W = 1 * pow(10,6); // Mass per Year (g/yr)
const double k = 0.25; // Volume per Mass per Year (m^0.5/g^0.5/yr)

// Steady-State Concentration of Pollutant
double pSS(double c,double t)
{
	double p = -c + t*(W/V) - t*((Q*c)/V) - pow(c,0.5)*k*t;
	return p;
}

// Problem Function
double mainF(double xm,double ym)
{
	double Y = pSS(xm,ym);
	return Y;
}

void mainR(double i0,double i1,double i2,double i3,double i4,double i5,int i6)
{
	printf("%2.0f  %2.3f  %2.2f  %2.3f  %2.3f  %2.3f   %2i \n",i0,i1,i2,i3,i4,i5,i6);
}

// Test Function
double testF(double xt)
{
	double G = pow(xt,2) - 5*xt + 6;
	return G;
}

void testR(double j0,double j1,double j2,double j3,double j4,int j5)
{
	printf("%2.1f %2.3f %2.3f %2.3f %2.3f %2i \n",j0,j1,j2,j3,j4,j5);
}

// Modified Secant Method
void ModSecant(int F,double C,double X0,double D,double Es,int Imax,int Ilim)
{
	int iter = 0; // Iteration
	double ea = Es + 1; // Absolute Error, actual (initial)
	double xr = X0; // Root of Function (initial)
	double xrold = 0; // Root of Function (previous)
	double fr = 0; // Function Value at Root (initial)

	do
	{
		xrold = xr;
		if (F == 0) xr = xrold - ((D*xrold*mainF(xrold,C))/(mainF(xrold+D*xrold,C)-mainF(xrold,C)));
		else if (F == 1) xr = xrold - ((D*xrold*testF(xrold))/(testF(xrold+D*xrold)-testF(xrold)));
		iter++;
		if (xr != 0) ea = abs((xr-xrold)/xr)*100;
		if (F == 0 && iter <= Ilim)
		{
			fr = mainF(xr,C);
			printf(" %2i   %2.1f  %2.3f %2.2f %2.3f  %2.3f  %2.3f \n",iter,C,D,X0,xr,fr,ea);
		}
	} while (ea >= Es && iter < Imax);
	if (ea < Es)
	{
		if (F == 0 && Ilim <= -1)
		{
			fr = mainF(xr,C);
			mainR(C,D,X0,xr,fr,ea,iter);
		}
		else if (F == 1)
		{
			fr = testF(xr);
			testR(X0,D,xr,fr,ea,iter);
		}
	}
	else if (iter >= Imax) std::cout << "No Solution due to Divergence" << std::endl;
}

int main()
{
	const double dT = 0.001; // Increment, Test
	const double dP = 0.5; // Increment, Problem
	const double es = 0.0005; // Error Criteria;
	const int imax = 100; // Maximum Number of Iterations

	std::cout << "Chapter 6 - Problem 6.18" << std::endl;
	std::cout << "============================================" << std::endl;
	std::cout << "Guess_initial, x0" << std::endl;
	std::cout << "Increment, d" << std::endl;
	std::cout << "Root of Function, xr" << std::endl;
	std::cout << "Function Value at Root, f(xr)" << std::endl;
	std::cout << "Error Criteria, es = " << es << std::endl;
	std::cout << "Actual Error, ea" << std::endl;
	std::cout << "Maximum Number of Iterations, imax = " << imax << std::endl;
	std::cout << "Iterations, iter" << std::endl;
	std::cout << "********************************************" << std::endl;

	std::cout << "Test Function" << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	std::cout << "x0  d     xr     f(xr)  ea    iter" << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	ModSecant(1,0,5,dT,es,imax,0);
	ModSecant(1,0,1,dT,es,imax,0);
	std::cout << "********************************************" << std::endl;

	std::cout << "Problem Function" << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	std::cout << "Volume, V(m^3) = " << V << std::endl;
	std::cout << "Volume per Year, Q(m^3/yr) = " << Q << std::endl;
	std::cout << "Mass per Year, W(g/yr) = " << W << std::endl;
	std::cout << "Volume per Mass per Year, k(m^0.5/g^0.5/yr) = " << k << std::endl;
	std::cout << "Time, t(yr)" << std::endl;
	std::cout << "cVariable_initial, c0(g/m^3)" << std::endl;
	std::cout << "cVaribale_root, cr(g/m^3)" << std::endl;
	std::cout << "Function Value at cVaribale_root, f(cr)" << std::endl;
	std::cout << "Iteration, iter." << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	std::cout << "iter. t     d     c0   cr     f(cr)  ea" << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	ModSecant(0,10,4,dP,es,imax,3);
	std::cout << "********************************************" << std::endl;
	std::cout << "t   d      c0    cr     f(cr)  ea     iter" << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	ModSecant(0,5,4,dP,es,imax,-1);
	ModSecant(0,10,4,dP,es,imax,-1);
	ModSecant(0,15,4,dP,es,imax,-1);
	ModSecant(0,20,4,dP,es,imax,-1);
	ModSecant(0,25,4,dP,es,imax,-1);
	ModSecant(0,30,4,dP,es,imax,-1);
	std::cout << "******************************************** \n" << std::endl;
	
	system("pause");
	return 0;
}