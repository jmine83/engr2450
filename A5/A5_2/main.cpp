// 03/25/2014 - ENGR 2450 - Meine, Joel
// Problems 21.3, 24.34

// Simpson's Rules (Equally-Spaced)

#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
using namespace std;

double Trap(double h,double f0,double f1)
{
	return((h/2) * (f0+f1));
}

double Simp38(double h,double f0,double f1,double f2,double f3)
{
	return(3*h * ((f0+3*(f1+f2)+f3)/8));
}

double Simp13m(double h,int n,vector<double> f)
{
	double sum = f[0];
	for (int i = 1; i < (n-2); i = i+2)
		sum = sum + 4*f[i] + 2*f[i+1];
	sum = sum + 4*f[n-1] + f[n];
	return(h * (sum/3));
}

double SimpInt(double h,double a,double b,int n,vector<double> f)
{
	double sum = 0;
	if (n == 1)
		sum = Trap(h,f[n-1],f[n]);
	else
	{
		int m = n;
		int odd = n%2;
		if (odd == 1 && n > 1)
		{
			sum = sum + Simp38(h,f[n-3],f[n-2],f[n-1],f[n]);
			m = n - 3;
		}
		if (m > 1)
			sum = sum + Simp13m(h,m,f);

	}
	return(sum);
}

double F(double x)
{
	return(1 - x - 4*pow(x,3) + 2*pow(x,5));
}

void funcEval(double a,double b,int n,double& h,double& I)
{
	h = (b-a)/n;
	vector<double> y;
	for (int i = 0; i <= n; i++)
	{
		y.push_back(F(a+i*h));
	}
	I = SimpInt(h,a,b,n,y);
	return;
}

void dataEval(double a,double b,int n,vector<double> y,double& I)
{
	double h = (b-a)/n;
	I = SimpInt(h,a,b,n,y);
	return;
}

int main()
{
	double h = 0; // Size of Subinterval, h

	// Problem 21.3 - Integral of a Function
	const double xi = -2; // Initial Value, xi
	const double xf = 4; // Final Value, xf
	double s = 0; // Result of Integration, s

	std::cout << "Chapter 21 - Problem 21.3" << std::endl;
	std::cout << "============================================" << std::endl;
	std::cout << "Initial Value, xi = " << setprecision(0) << xi << std::endl;
	std::cout << "Final Value, xf = " << setprecision(0) << xf << std::endl;
	std::cout << "Number of Subintervals, n" << std::endl;
	std::cout << "Size of Subinterval, h" << std::endl;
	std::cout << "Result of Integration, s" << std::endl;
	std::cout << "********************************************" << std::endl;
	std::cout << " n   h    s" << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	const int n_size = 5;
	int n[] = {2,3,4,10,20};
	for (int i = 0; i < n_size; i++)
	{
		funcEval(xi,xf,n[i],h,s);
		cout << setw(3) << n[i];
		cout << setw(5) << setprecision(3) << h;
		cout << setw(11) << setprecision(9) << s << endl;
	}
	std::cout << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	cout << "\n";

	// Problem 24.34 - Integral for Tabulated Data
	vector<double> t = {0.0,0.2,0.4,0.6,0.8,1.0,1.2}; // Time (s), t
	vector<double> I = {0.2,0.3683,0.3819,0.2282,0.0486,0.0082,0.1441}; // Current (mA), I(t)
	double ti = t.front(); // Initial Time (s), ti
	double tf = t.back(); // Final Time (s), tf
	const double C = 10 * pow(10,-6); // Capacitance (F), C
	double Is = 0; // Total Current (mA), Is

	std::cout << "Chapter 24 - Problem 24.34" << std::endl;
	std::cout << "============================================" << std::endl;
	std::cout << "Initial Time (s), ti = " << setprecision(0) << ti << std::endl;
	std::cout << "Final Time (s), tf = " << setprecision(0) << tf << std::endl;
	std::cout << "Voltage (V), V = 1/C * Is" << std::endl;
	std::cout << "Capacitance (F), C = " << setprecision(7) << C << std::endl;
	std::cout << "Total Current (mA), Is" << std::endl;
	std::cout << "********************************************" << std::endl;

	dataEval(ti,tf,I.size()-1,I,Is);
	double V = (1/C) * Is; // Voltage (V), V
	std::cout << "Voltage (V), V = " << setprecision(7) << V << std::endl;
	std::cout << "Total Current (mA), Is = " << setprecision(7) << Is << std::endl;

	std::cout << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	cout << "\n";

	system("pause");
	return 0;
}