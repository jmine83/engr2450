// 03/05/2014 - ENGR 2450 - Meine, Joel
// Problems 18.6, 20.36

// Newton's Polynomial Interpolation

#include <iostream>
#include <iomanip>
#include <math.h>
using namespace std;

const int n = 4; // Polynomial Order, max

void NewInt(double x[],double y[],double xi,double yint[],double ea[])
{
	double yint2 = 0;
	double fdd[n+1][n+1];
	for (int i = 0; i <= n; i++)
		fdd[i][0] = y[i];
	for (int j = 1; j <= n; j++)
	{
		for (int i = 0; i <= n - j; i++)
			fdd[i][j] = (fdd[i+1][j-1] - fdd[i][j-1])/(x[i+j] - x[i]);
	}
	double xterm = 1;
	yint[0] = fdd[0][0];
	for (int order = 1; order <= n; order++)
	{
		xterm = xterm * (xi - x[order-1]);
		yint2 = yint[order-1] + fdd[0][order] * xterm;
		ea[order-1] = yint2 - yint[order-1];
		yint[order] = yint2;
	}
	ea[n] = 0;
}

int main()
{
	// Problem 18.6 - Newton's Polynomial Interpolation, Order 1 to 4
	double x[n+1] = {5,3,7,2,8};
	double y[n+1] = {99,19,291,6,444};
	double xint = 4; // Value of Interest
	double yint[n+1]; // Solution of Interest
	double ea1[n+1];
	NewInt(x,y,xint,yint,ea1);

	std::cout << "Chapter 18 - Problem 18.6" << std::endl;
	std::cout << "============================================" << std::endl;
	std::cout << "Value of Interest, xint = " << xint << std::endl;
	std::cout << "Solution of Interest, yint" << std::endl;
	std::cout << "Actual Error, ea" << std::endl;
	std::cout << "********************************************" << std::endl;
	std::cout << " x     y" << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	for (int i = 0; i < n + 1; i++)
	{
		cout << setw(2) << x[i];
		cout << setw(6) << y[i] << endl;
	}
	std::cout << "********************************************" << std::endl;
	std::cout << " yint   ea" << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	for (int i = 0; i < n + 1; i++)
	{
		cout << setw(5) << setprecision(5) << yint[i];
		cout << setw(5) << setprecision(5) << ea1[i] << endl;
	}
	std::cout << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	cout << "\n";

	// Problem 20.36 - Newton's Polynomial Interpolation, Order 1 to 4
	double I[n+1] = {1.25,0.75,1.5,0.25,2.0}; // Current (A), I
	double V[n+1] = {0.70,-0.6,1.88,-0.45,6.0}; // Voltage (V), V
	double Iint = 1.15; // Value of Interest
	double Vint[n+1]; // Solution of Interest
	double ea2[n+1];
	NewInt(I,V,Iint,Vint,ea2);

	std::cout << "Chapter 20 - Problem 20.36" << std::endl;
	std::cout << "============================================" << std::endl;
	std::cout << "Current (A), I" << std::endl;
	std::cout << "Voltage (V), V" << std::endl;
	std::cout << "Value of Interest, Iint = " << Iint << std::endl;
	std::cout << "Solution of Interest, Vint" << std::endl;
	std::cout << "Actual Error, ea" << std::endl;
	std::cout << "********************************************" << std::endl;
	std::cout << " I      V" << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	for (int i = 0; i < n + 1; i++)
	{
		cout << setw(2) << I[i];
		cout << setw(8) << V[i] << endl;
	}
	std::cout << "********************************************" << std::endl;
	std::cout << " Vint       ea" << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	for (int i = 0; i < n + 1; i++)
	{
		cout << setw(9) << setprecision(6) << Vint[i];
		cout << setw(14) << setprecision(6) << ea2[i] << endl;
	}
	std::cout << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

	cout << "\n";
	system("pause");
	return 0;
}