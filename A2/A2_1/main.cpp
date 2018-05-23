// 02/05/2014 - ENGR 2450 - Meine, Joel
// Chapter 5 - Problem 5.18

#include <iostream>
#include <math.h>
using namespace std;

// Oxygen Saturation in Freshwater, osf(mg/L)
const double Tl = 0; // Temperature_lower, Tl(C)
const double Tu = 40; // Temperature_upper, Tu(C)

double Osf(double T,double osf)
{
	double Ta = T + 273.15; // Absolute Temperature, Ta(K); Temperature, T(C)
	double O = -(log(osf)) - 139.34411 + ((1.575701*pow(10,5)) / (pow(Ta,1))) - ((6.642308*pow(10,7)) / (pow(Ta,2))) + ((1.243800*pow(10,10)) / (pow(Ta,3))) - ((8.621949*pow(10, 11)) / (pow(Ta,4)));
	return O;
}

// Problem Function
double mainF(double xm,double ym)
{
	double Y = Osf(xm,ym);
	return Y;
}

void mainR(double i0,double i1,double i2,double i3,double i4,double i5,int i6)
{
	printf("%2.0f  %2.0f   %2.0f  %2.4f  %2.4f  %2.4f  %2i \n",i0,i1,i2,i3,i4,i5,i6);
}

// Test Function
double testF(double xt)
{
	double G = pow(xt,2) - 5*xt + 6;
	return G;
}

void testR(double j0,double j1,double j2,double j3,double j4,int j5)
{
	printf("%2.1f  %2.1f  %2.4f  %2.4f  %2.4f  %2i \n",j0,j1,j2,j3,j4,j5);
}

// Bisection Method
void BiSect(int F,double C,double xl,double xu,double Ead,int Imax)
{
	double xli = 0; // Bracket_lower (initial)
	double xui = 0; // Bracket_upper (initial)
	int iter = 0; // Iteration
	double ea = Ead + 1; // Absolute Error, actual (initial)
	double xr = 0; // Root of Function (initial)
	double xrold = 0; // Root of Function (previous)

	double fl = 0; // Function Value at Bracket_lower (initial)
	if (F == 0) fl = mainF(xl,C);  // Function Value at Bracket_lower (actual function)
	else if (F == 1) fl = testF(xl);  // Function Value at Bracket_lower (test function)
	double fu = 0; // Function Value at Bracket_upper (initial)
	if (F == 0) fu = mainF(xu,C);  // Function Value at Bracket_upper (actual function)
	else if (F == 1) fu = testF(xu);  // Function Value at Bracket_upper (test function)
	double fr = 0; // Function Value at Root (initial)

	do
	{
		if (iter == 0)
		{
			xli = xl;
			xui = xu;
		}
		xrold = xr;
		xr = (xl + xu) / 2;
		ea = abs(xr-xrold);
		if (F == 0) fr = mainF(xr,C);
		else if (F == 1) fr = testF(xr);
		iter++;
		double test = fl * fr;
		if (test < 0) xu = xr;
		else if (test > 0)
		{
			xl = xr;
			fl = fr;
		}
		else ea = 0;
	} while (ea >= Ead && iter < Imax && fl*fu < 0);
	if (ea < Ead)
	{
		if (F == 0)
		{
			fr = mainF(xr,C);
			mainR(C,xli,xui,xr,fr,ea,iter);
		}
		else if (F == 1)
		{
			fr = testF(xr);
			testR(xli,xui,xr,fr,ea,iter);
		}
	}
	else if (iter >= Imax) std::cout << "No Solution due to Divergence" << std::endl;
	else if (fl*fu >= 0) std::cout << "Invalid Values for xl and xu" << std::endl;
}

int main()
{
	const double ead = 0.05; // Absolute Error
	const int imax = 100; // Maximum Number of Iterations
	
	std::cout << "Chapter 5 - Problem 5.18" << std::endl;
	std::cout << "============================================" << std::endl;
	std::cout << "Bracket_lower, xl" << std::endl;
	std::cout << "Bracket_upper, xu" << std::endl;
	std::cout << "Root of Function, xr" << std::endl;
	std::cout << "Function Value at Root, f(xr)" << std::endl;
	std::cout << "Absolute Error, ead = " << ead << std::endl;
	std::cout << "Actual Error, ea" << std::endl;
	std::cout << "Maximum Number of Iterations, imax = " << imax << std::endl;
	std::cout << "Iterations, iter" << std::endl;
	std::cout << "********************************************" << std::endl;

	std::cout << "Test Function" << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	std::cout << "xl   xu   xr      f(xr)   ea      iter" << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	BiSect(1,0,1.0,2.5,ead,imax);
	BiSect(1,0,2.5,4.0,ead,imax);
	std::cout << "********************************************" << std::endl;

	std::cout << "Problem Function" << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	std::cout << "Oxygen Saturation in Freshwater, osf(mg/L)" << std::endl;
	std::cout << "Temperature_lower, Tl(C)" << std::endl;
	std::cout << "Temperature_upper, Tu(C)" << std::endl;
	std::cout << "Temperature_root, Tr(C)" << std::endl;
	std::cout << "Function Value at Temperature_root, f(Tr)" << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	std::cout << "osf  Tl  Tu  Tr       f(Tr)   ea     iter" << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	BiSect(0,8,Tl,Tu,ead,imax); // osf = 8
	BiSect(0,10,Tl,Tu,ead,imax); // osf = 10
	BiSect(0,12,Tl,Tu,ead,imax); // osf = 12
	std::cout << "******************************************** \n" << std::endl;

	system("pause");
	return 0;
}