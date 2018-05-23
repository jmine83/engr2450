// 02/05/2014 - ENGR 2450 - Meine, Joel
// Chapter 5 - Problem 5.22

#include <iostream>
#include <math.h>
using namespace std;

// Population of Urban and Suburban Areas
const double Pumax = 75000; // Population, Urban (maximum)
const double ku = 0.045; // Population Rate of Change, Urban
const double Pumin = 100000; // Population, Urban (minimum)
const double Psmax = 300000; // Population, Suburban (maximum)
const double ks = 0.08; // Population Rate of Change, Suburban
const double P0 = 10000; // Population (initial)

const double Tl = 0; // Time_lower, Tl(yr)
const double Tu = 100; // Time_upper, TU(yr)

double P(double T,double Pc)
{
	double Pu = Pumax * exp(-ku*T) + Pumin; // Population, Urban
	double Ps = Psmax / (1 + (Psmax/(P0-1)) * exp(-ks*T)); // Population, Suburban
	double Pt = Pc - (Ps/Pu); // Suburb-to-Urban Population Comparison
	return Pt;	
}

// Problem Function
double mainF(double xm,double ym)
{
	double Y = P(xm,ym);
	return Y;
}

void mainR(double i0,double i1,double i2,double i3,double i4,double i5,int i6)
{
	printf("%2.1f  %2.0f   %2.0f  %2.4f  %2.4f  %2.4f  %2i \n",i0,i1,i2,i3,i4,i5,i6);
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

// Modified False-Position Method
void ModFalsePos(int F,double C,double xl,double xu,double Es,int Imax)
{
	double xli = 0; // Bracket_lower (initial)
	double xui = 0; // Bracket_upper (initial)
	int iter = 0; // Iteration
	int il = 0; // Iteration_lower
	int iu = 0; // Iteration_upper
	double ea = Es + 1; // Absolute Error, actual (initial)
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
		xr = xu - fu * (xl-xu)/(fl-fu);
		if (F == 0) fr = mainF(xr,C);
		else if (F == 1) fr = testF(xr);
		iter++;
		if (xr != 0) ea = abs((xr-xrold)/xr)*100;
		double test = fl * fr;
		if (test < 0)
		{
			xu = xr;
			if (F == 0) fu = mainF(xu,C);
			else if (F == 1) fu = testF(xu);
			iu = 0;
			il++;
			if (il >= 2) fl = fl/2;
		}
		else if (test > 0)
		{
			xl = xr;
			if (F == 0) fl = mainF(xl,C);
			else if (F == 1) fl = testF(xl);
			il = 0;
			iu++;
			if (iu >= 2) fu = fu/2;
		}
		else ea = 0;
	} while (ea >= Es && iter < Imax && fl*fu < 0);
	if (ea < Es)
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
	const double es = 0.0005; // Error Criteria
	const int imax = 100; // Maximum Number of Iterations

	std::cout << "Chapter 5 - Problem 5.22" << std::endl;
	std::cout << "============================================" << std::endl;
	std::cout << "Bracket_lower, xl" << std::endl;
	std::cout << "Bracket_upper, xu" << std::endl;
	std::cout << "Root of Function, xr" << std::endl;
	std::cout << "Function Value at Root, f(xr)" << std::endl;
	std::cout << "Error Criteria, es = " << es << std::endl;
	std::cout << "Actual Error, ea" << std::endl;
	std::cout << "Maximum Number of Iterations, imax = " << imax << std::endl;
	std::cout << "Iterations, iter" << std::endl;
	std::cout << "********************************************" << std::endl;

	std::cout << "Test Function" << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	std::cout << "xl   xu   xr      f(xr)   ea      iter" << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	ModFalsePos(1,0,1.0,2.5,es,imax);
	ModFalsePos(1,0,2.5,4.0,es,imax);
	std::cout << "********************************************" << std::endl;

	std::cout << "Problem Function" << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	std::cout << "Population, Urban (maximum), Pumax = " << Pumax << std::endl;
	std::cout << "Population Rate of Change (Urban), ku = " << ku << std::endl;
	std::cout << "Population, Urban (minimum), Pumin = " << Pumin << std::endl;
	std::cout << "Population, Suburban (maximum), Psmax = " << Psmax << std::endl;
	std::cout << "Population Rate of Change (Suburban), ks = " << ks << std::endl;
	std::cout << "Population (initial), P0 = " << P0 << std::endl;
	std::cout << "Suburb-to-Urban Population Comparison, Pc" << std::endl;
	std::cout << "Time_lower, Tl(yr)" << std::endl;
	std::cout << "Time_upper, Tu(yr)" << std::endl;
	std::cout << "Time_root, Tr(yr)" << std::endl;
	std::cout << "Function Value at Time_root, f(Tr)" << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	std::cout << "Pc    Tl  Tu  Tr        f(Tr)   ea     iter" << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	ModFalsePos(0,1.2,Tl,Tu,es,imax); // Pc = 1.2; Pc > 1.% (% larger than) && Pc < 1.% (% smaller than)
	std::cout << "******************************************** \n" << std::endl;

	system("pause");
	return 0;
}