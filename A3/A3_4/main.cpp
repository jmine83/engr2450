// 02/19/2014 - ENGR 2450 - Meine, Joel
// Chapter 10 - Problem 10.19

// Inverse Matrix

#include <iostream>
#include <iomanip>
#include <math.h>
using namespace std;

const int n = 3;
const double tol = 0.0001;

void Pivot(double a[][n],int o[],double s[],int k)
{
	int p = k;
	double dummy = 0;
	double big = abs(a[o[k]][k] / s[o[k]]);
	for (int ii = k+1; ii < n; ii++)
	{
		dummy = abs(a[o[ii]][k] / s[o[ii]]);
		if (dummy > big)
		{
			big = dummy;
			p = ii;
		}
	}
	dummy = o[p];
	o[p] = o[k];
	o[k] = dummy;
}

void Decompose(double a[][n],double s[],int o[],int er)
{
	for (int i = 0; i < n; i++)
	{
		o[i] = i;
		s[i] = abs(a[i][0]);
		for (int j = 1; j < n; j++)
		{
			if (abs(a[i][j]) > s[i]) s[i] = abs(a[i][j]);
		}
	}
	for (int k = 0; k < (n-1); k++)
	{
		Pivot(a,o,s,k);
		if (abs(a[o[k]][k] / s[o[k]]) < tol)
		{
			er = -1;
			break;
		}
		for (int i = k+1; i < n; i++)
		{
			double factor = a[o[i]][k] / a[o[k]][k];
			a[o[i]][k] = factor;
			for (int j = k+1; j < n; j++) a[o[i]][j] = a[o[i]][j] - factor * a[o[k]][j];
		}
	}
	if (abs(a[o[n-1]][n-1] / s[o[n-1]]) < tol) er = -1;
}

void Substitute(double a[][n],int o[],double b[],double x[])
{
	for (int i = 1; i < n; i++)
	{
		double sum = b[o[i]];
		for (int j = 0; j <= (i-1); j++) sum = sum - a[o[i]][j] * b[o[j]];
		b[o[i]] = sum;
	}
	x[n-1] = b[o[n-1]] / a[o[n-1]][n-1];
	for (int i = n-2; i >= 0; i--)
	{
		double sum = 0;
		for (int j = i+1; j < n; j++) sum = sum + a[o[i]][j] * x[j];
		x[i] = (b[o[i]] - sum) / a[o[i]][i];
	}
}

void MatrixInverse(double a[][n],double ai[][n],int er)
{
	int o[n];
	double s[n];
	double b[n];
	double x[n];
	er = 0;
	Decompose(a,s,o,er);
	if (er == 0)
	{
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				if (i == j) b[j] = 1;
				else b[j] = 0;
			}
			Substitute(a,o,b,x);
			for (int j = 0; j < n; j++) ai[j][i] = x[j];
		}
	}
	else std::cout << "Ill-Conditioned System" << std::endl;
}

void copymatrix(double A[][n],double B[][n],int n,int m)
{
	int i; int j;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
		{
			B[i][j] = A[i][j];
		}
	}
}

void multiply_matrices(double A[][n],double B[][n],double C[][n],int n,int m,int p)
{
	int i; int j; int k;
	for (i = 0; i < n; i++){
		for (j = 0; j < m; j++){
			C[i][j] = 0;
			for (k = 0; k < p; k++){
				C[i][j] = C[i][j] + A[i][k] * B[k][j];
			};
		};
	};
}

void printmatrix(double A[][n],int n,int m,int p)
{
	int i; int j;
	for (i = 0; i < n; i++){
		for (j = 0; j < m; j++){
			cout << setw(10) << setprecision(p) << A[i][j] << " ";
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
	cout << "A = \n"; printmatrix(A,n,n,0); // Matrix A, Print
	double AA[n][n]; // Matrix AA, Initialize
	copymatrix(A,AA,n,n); // AA = A
	cout << "AA = \n"; printmatrix(AA,n,n,0); // Matrix AA, Print
	double AI[n][n]; // Matrix AI, Intialize
	MatrixInverse(AA,AI,0); // AI = AA^-1
	cout << "AI = \n"; printmatrix(AI,n,n,3); // Matrix AI, Print
	double b1[n][n]; // Matrix b1, Initialize
	multiply_matrices(A,AI,b1,n,n,n); // b1 = A * AI
	int b2[n][n]; // Matrix b2, Initialize
	double B[n][n]; // Matrix B, Initialize
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			b2[i][j] = b1[i][j];
			B[i][j] = b2[i][j];
		}
	}
	cout << "B = \n"; printmatrix(B,n,n,0); // Matrix B, Print
	
	cout << "\n";
	system("pause");
	return 0;
}