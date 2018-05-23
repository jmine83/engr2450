// 03/25/2014 - ENGR 2450 - Meine, Joel
// Problems 21.22, 24.4

// Trapezoidal Rule & Simpson's Rules (Unequally-Spaced)

#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
using namespace std;

double Trapun(vector<double> x,vector<double> y,int n)
{
	double sum = 0;
	for (int i = 1; i <= n; i++)
		sum = sum + (((x[i]-x[i-1])*(y[i-1]+y[i]))/2);
	return(sum);
}

double Trap(double h,double f0,double f1)
{
	return(h * ((f0+f1)/2));
}

double Simp38(double h,double f0,double f1,double f2,double f3)
{
	return(3*h * ((f0+3*(f1+f2)+f3)/8));
}

double Simp13(double h,double f0,double f1,double f2)
{
	return(2*h * ((f0+4*f1+f2)/6));
}

double Uneven(int n,vector<double> x,vector<double> f)
{
	double h = x[1] - x[0];
	int k = 1;
	double sum = 0;
	double hf = 0;
	for (int j = 1; j <= n; j++)
	{
		if (j == n)
			hf = x[0] - x[j];
		else
			hf = x[j+1] - x[j];
		if (abs(h-hf) < .000001)
		{
			if (k == 3)
			{
				sum = sum + Simp13(h,f[j-3],f[j-2],f[j-1]);
				k = k - 1;
			}
			else
				k = k + 1;
		}
		else
		{
			if (k == 1)
			{
				sum = sum + Trap(h,f[j-1],f[j]);
			}
			else
			{
				if (k == 2)
				{
					sum = sum + Simp13(h,f[j-2],f[j-1],f[j]);
				}
				else
				{
					sum = sum + Simp38(h,f[j-3],f[j-2],f[j-1],f[j]);
				}
				k = 1;
			}
		}
		h = hf;
	}
	return(sum);
}

double F(double x)
{
	return(1 - x - 4*pow(x,3) + 2*pow(x,5));
}

int main()
{
	// Problem 21.22 - Trapezoidal Rule (Unequally-Spaced)
	vector<double> V = {0.5,2,3,4,6,8,10,11}; // Volume (m^3), V
	vector<double> p = {336,294.4,266.4,260.8,260.5,249.6,193.6,165.6}; // Pressure (kPa), p(V)
	double ps = 0; // Total Pressure (kPa), ps
	int n = V.size()-1;

	std::cout << "Chapter 21 - Problem 21.22" << std::endl;
	std::cout << "============================================" << std::endl;
	std::cout << "Work (kJ), W = ps" << std::endl;
	std::cout << "Total Pressure (kPa), ps" << std::endl;
	std::cout << "Number of Subintervals, n = " << n << std::endl;
	std::cout << "********************************************" << std::endl;

	double W = Trapun(V,p,n);
	std::cout << "Work (kJ), W = " << setprecision(7) << W << std::endl;

	std::cout << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	cout << "\n";
	
	// Problem 24.4 - Simpson's Rules (Unequally-Spaced)
	vector<double> t = {0,10,20,30,35,40,45,50}; // Time (min), t
	vector<double> c = {10,35,55,52,40,37,32,34}; // Mass Concentration (mg/m^3), c
	const double Q = 4; // Flow Rate Constant (m^3/min), Q
	vector<double> Qc;
	for (int i = 0; i < c.size(); i++)
		Qc.push_back(Q*c[i]);
	double cs = 0; // Total Mass Concentration (mg/m^3), cs
	int m = t.size()-1;

	std::cout << "Chapter 24 - Problem 24.4" << std::endl;
	std::cout << "============================================" << std::endl;
	std::cout << "Mass (mg), M = Q * cs" << std::endl;
	std::cout << "Flow Rate Constant (m^3/min), Q = " << Q << std::endl;
	std::cout << "Total Mass Concentration (mg/m^3), cs" << std::endl;
	std::cout << "Number of Subintervals, n = " << m << std::endl;
	std::cout << "********************************************" << std::endl;

	double M = Uneven(m,t,Qc);
	std::cout << "Mass (mg), M = " << setprecision(7) << M << std::endl;

	std::cout << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	cout << "\n";

	system("pause");
	return 0;
}