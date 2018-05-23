// 02/19/2014 - ENGR 2450 - Meine, Joel
// Chapter 11 - Problem 11.26

// Gauss-Seidel with Relaxation

#include <iostream>
#include <iomanip>
#include <math.h>
using namespace std;

const double es = 0.05;
const int imax = 50;
const int n = 3;

void Gseid(double a[][n],double b[],double x[],int imax,double es,double lambda)
{
	double sum = 0;
	for (int i = 0; i < n; i++)
	{
		double dummy = a[i][i];
		for (int j = 0; j < n; j++) a[i][j] = a[i][j] / dummy;
		b[i] = b[i] / dummy;
	}
	for (int i = 0; i < n; i++)
	{
		sum = b[i];
		for (int j = 0; j < n; j++)
		{
			if (i != j) sum = sum - a[i][j] * x[j];
		}
		x[i] = sum;
	}
	int iter = 1;
	int sentinel = 0;
	do
	{
		sentinel = 1;
		for (int i = 0; i < n; i++)
		{
			double old = x[i];
			sum = b[i];
			for (int j = 0; j < n; j++)
			{
				if (i != j) sum = sum - a[i][j] * x[j];
			}
			x[i] = lambda*sum + (1 - lambda)*old;
			if (sentinel == 1 && x[i] != 0)
			{
				double ea = abs((x[i]-old)/x[i])*100;
				if (ea > es) sentinel = 0;
			}
		}
		iter++;
	} while (sentinel != 1 || iter < imax);
	if (sentinel == 1) ;
	else if (iter >= imax) std::cout << "No Solution due to Divergence" << std::endl;
}

void printmatrix(double A[][n],int n,int m)
{
	int i; int j;
	for (i = 0; i < n; i++){
		for (j = 0; j < m; j++){
			cout << setw(10) << A[i][j] << " ";
		}; cout << endl;
	};
}

void printvector(double u[],int n)
{
	int i;
	for (i = 0; i < n; i++){
		cout << setw(10) << u[i]; cout << endl;
	};
}

int main()
{
	const double L1 = 0.5;
	const double L2 = 1.5;

	double A[n][n] = {{3,-0.1,-0.2},{0.1,7,-0.3},{0.3,-0.2,10}}; // Matrix A, Initialize
	cout << "A = \n"; printmatrix(A,n,n); // Matrix A, Print
	double b[n] = {7.85,-19.3,71.4}; // Vector b, Intialize
	cout << "b = \n"; printvector(b,n); // Vector b, Print
	double X1[n]; // Vector X1, Intialize
	Gseid(A,b,X1,imax,es,L1); // lambda = 0.5
	std::cout << "\nlambda = " << L1 << std::endl;
	cout << "X = \n"; printvector(X1,n); // Vector X1, Print
	double X2[n]; // Vector X2, Initialize
	Gseid(A,b,X2,imax,es,L2); // lambda = 1.5
	std::cout << "\nlambda = " << L2 << std::endl;
	cout << "X = \n"; printvector(X2,n); // Vector X2, Print
	
	cout << "\n";
	system("pause");
	return 0;
}