// 02/19/2014 - ENGR 2450 - Meine, Joel
// Chapter 9 - Problem 9.18

// Gauss Elimination with Partial Pivoting

#include <iostream>
#include <iomanip>
#include <math.h>
using namespace std;

const int n = 3;
const double tol = 0.0001;

void Pivot(double a[][n],double b[],double s[],int k)
{
	int p = k;
	double dummy = 0;
	double big = abs(a[k][k] / s[k]);
	for (int ii = k+1; ii < n; ii++)
	{
		dummy = abs(a[ii][k] / s[ii]);
		if (dummy > big)
		{
			big = dummy;
			p = ii;
		}
	}
	if (p != k)
	{
		for (int jj = k; jj < n; jj++)
		{
			dummy = a[p][jj];
			a[p][jj] = a[k][jj];
			a[k][jj] = dummy;
		}
		dummy = b[p];
		b[p] = b[k];
		b[k] = dummy;
		dummy = s[p];
		s[p] = s[k];
		s[k] = dummy;
	}
}

void Eliminate(double a[][n],double s[],double b[],int er)
{
	for (int k = 0; k < (n-1); k++)
	{
		Pivot(a,b,s,k);
		if (abs(a[k][k] / s[k]) < tol)
		{
			er = -1;
			break;
		}
		for (int i = k+1; i < n; i++)
		{
			double factor = a[i][k] / a[k][k];
			for (int j = k+1; j < n; j++) a[i][j] = a[i][j] - factor * a[k][j];
			b[i] = b[i] - factor * b[k];
		}
	}
	if (abs(a[n-1][n-1] / s[n-1]) < tol) er = -1;
}

void Substitute(double a[][n],double b[],double x[])
{
	x[n-1] = b[n-1] / a[n-1][n-1];
	for (int i = n-2; i >= 0; i--)
	{
		double sum = 0;
		for (int j = i+1; j < n; j++) sum = sum + a[i][j] * x[j];
		x[i] = (b[i] - sum) / a[i][i];
	}
}

void Gauss(double a[][n],double b[],double x[],int er)
{
	double s[n];
	er = 0;
	for (int i = 0; i < n; i++)
	{
		s[i] = abs(a[i][0]);
		for (int j = 1; j < n; j++)
		{
			if (abs(a[i][j]) > s[i]) s[i] = abs(a[i][j]);
		}
	}
	Eliminate(a,s,b,er);
	if (er != -1) Substitute(a,b,x);
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
	double A[n][n] = {{1,2,-1},{5,2,2},{-3,5,-1}}; // Matrix A, Initialize
	cout << "A = \n"; printmatrix(A,n,n); // Matrix A, Print
	double b[n] = {2,9,1}; // Vector b, Intialize
	cout << "b = \n"; printvector(b,n); // Vector b, Print
	double x[n]; // Vector x, Intialize
	Gauss(A,b,x,0);
	cout << "x = \n"; printvector(x,n); // Vector x, Print
	
	cout << "\n";
	system("pause");
	return 0;
}