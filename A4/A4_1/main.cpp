// 03/05/2014 - ENGR 2450 - Meine, Joel
// Problems 17.4, 17.7

// Linear Regression

#include <iostream>
#include <iomanip>
#include <math.h>
using namespace std;

void Regress(double x[],double y[],int n,double& a1,double& a0,double& syx,double& r2,double& r)
{
	double sumx = 0; double sumxy = 0; double st = 0;
	double sumy = 0; double sumx2 = 0; double sr = 0;
	for (int i = 0; i < n; i++)
	{
		sumx = sumx + x[i];
		sumy = sumy + y[i];
		sumxy = sumxy + x[i]*y[i];
		sumx2 = sumx2 + x[i]*x[i];
	}
	double xm = sumx/n;
	double ym = sumy/n;
	a1 = (n*sumxy - sumx*sumy)/(n*sumx2 - sumx*sumx);
	a0 = ym - a1*xm;
	for (int i = 0; i < n; i++)
	{
		st = st + pow(y[i] - ym,2);
		sr = sr + pow(y[i] - a1*x[i] - a0,2);
	}
	syx = sqrt(sr/(n - 2));
	r2 = (st - sr)/st;
	r = sqrt(r2);
}

void FitData(double x[],int n,double a1,double a0,double yf[],int F)
{
	double alpha = 0; double beta = 0;
	for (int i = 0; i < n; i++)
	{
		if (F == 1) // Linear Equation
		{
			yf[i] = a0 + a1*x[i];
		}
		else if (F == 2) // Saturation-Growth-Rate Equation
		{
			alpha = 1/a0; beta = alpha*a1;
			yf[i] = alpha*(x[i]/(beta + x[i]));
		}
		else if (F == 3) // Power Equation
		{
			alpha = pow(10,a0); beta = a1;
			yf[i] = alpha*pow(x[i],beta);
		}
	}
}

int main()
{
	// Problem 17.4 - Least-Squares Regression
	const int n = 11;
	const int m = 11;
	double x1[n] = {6,7,11,15,17,21,23,29,29,37,39};
	double y1[m] = {29,21,29,14,21,15,7,7,13,0,3};
	double yf1[n];

	// Problem 17.7 - Data Fitting Types
	const int k = 7;
	const int j = 7;
	double x2[n] = {0.75,2,3,4,6,8,8.5};
	double y2[m] = {1.2,1.95,2,2.4,2.4,2.7,2.6};
		
	double yfA[n];
	double yfB[n];

	if (n == m || k == j)
	{
		// Problem 17.4 - Linear Equation
		double a1 = 0; // Slope
		double a0 = 0; // Intercept
		double syx = 0; // Standard Error of the Estimate
		double r2 = 0; // Coefficient of Determination
		double r = 0; // Correlation Coefficient
		
		Regress(x1,y1,n,a1,a0,syx,r2,r);

		std::cout << "Chapter 17 - Problem 17.4" << std::endl;
		std::cout << "============================================" << std::endl;
		std::cout << "Number of Data Points, n = " << n << std::endl;
		std::cout << "Slope, a1 = " << setprecision(4) << a1 << std::endl;
		std::cout << "Intercept, a0 = " << setprecision(4) << a0 << std::endl;
		std::cout << "Standard Error of Estimate, syx = " << setprecision(4) << syx << std::endl;
		std::cout << "Coefficient of Determination, r2 = " << setprecision(4) << r2 << std::endl;
		std::cout << "Correlation Coefficient, r = " << setprecision(4) << r << std::endl;
		std::cout << "********************************************" << std::endl;

		FitData(x1,n,a1,a0,yf1,1); // Linear Equation

		std::cout << "Linear Equation" << std::endl;
		std::cout << "--------------------------------------------" << std::endl;
		std::cout << " x   y    yf" << std::endl;
		std::cout << "--------------------------------------------" << std::endl;
		for (int i = 0; i < n; i++)
		{
			cout << setw(3) << x1[i];
			cout << setw(4) << y1[i];
			cout << setw(9) << setprecision(5) << yf1[i] << endl;
		}
		std::cout << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
		cout << "\n";
		
		// Problem 17.7 - Saturation-Growth-Rate Equation
		a1 = 0; // Slope
		a0 = 0; // Intercept
		syx = 0; // Standard Error of the Estimate
		r2 = 0; // Coefficient of Determination
		r = 0; // Correlation Coefficient
		
		Regress(x2,y2,k,a1,a0,syx,r2,r);

		std::cout << "Chapter 17 - Problem 17.7" << std::endl;
		std::cout << "============================================" << std::endl;
		std::cout << "Number of Data Points, n = " << k << std::endl;
		std::cout << "Slope, a1 = " << setprecision(4) << a1 << std::endl;
		std::cout << "Intercept, a0 = " << setprecision(4) << a0 << std::endl;
		std::cout << "Standard Error of Estimate, syx = " << setprecision(4) << syx << std::endl;
		std::cout << "Coefficient of Determination, r2 = " << setprecision(4) << r2 << std::endl;
		std::cout << "Correlation Coefficient, r = " << setprecision(4) << r << std::endl;
		std::cout << "********************************************" << std::endl;
		
		FitData(x2,k,a1,a0,yfA,2); // Saturation-Growth-Rate Equation
		
		std::cout << "Saturation-Growth-Rate Equation" << std::endl;
		std::cout << "--------------------------------------------" << std::endl;
		std::cout << " x     y    yf" << std::endl;
		std::cout << "--------------------------------------------" << std::endl;
		for (int i = 0; i < k; i++)
		{
			cout << setw(4) << x2[i];
			cout << setw(6) << y2[i];
			cout << setw(9) << setprecision(5) << yfA[i] << endl;
		}
		std::cout << "********************************************" << std::endl;

		FitData(x2,k,a1,a0,yfB,3); // Power Equation

		std::cout << "Power Equation" << std::endl;
		std::cout << "--------------------------------------------" << std::endl;
		std::cout << " x     y     yf" << std::endl;
		std::cout << "--------------------------------------------" << std::endl;
		for (int i = 0; i < k; i++)
		{
			cout << setw(4) << x2[i];
			cout << setw(6) << y2[i];
			cout << setw(9) << setprecision(5) << yfB[i] << endl;
		}
		std::cout << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	}
	else std::cout << "Regression is Not Possible" << std::endl;

	cout << "\n";
	system("pause");
	return 0;
}