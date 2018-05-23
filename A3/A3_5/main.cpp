// 02/19/2014 - ENGR 2450 - Meine, Joel
// Chapter 11 - Problem 11.25

// Cholesky Decomposition

#include <iostream>
#include <iomanip>
#include <math.h>
using namespace std;

const int n = 3;

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

void printmatrix(double A[][n],int n,int m, int s)
{
	int i; int j;
	for (i = 0; i < n; i++){
		for (j = 0; j < m; j++){
			if (i == 0 && j == 1 && s == 1) A[i][j] = 0;
			else if (i == 0 && j == 2 && s == 1) A[i][j] = 0;
			else if (i == 1 && j == 2 && s == 1) A[i][j] = 0;
			cout << setw(10) << A[i][j] << " ";
		}; cout << endl;
	};
}

void cLudecomp(double a[][n],double b[][n])
{
	double sum = 0;
	for (int k = 0; k < n; k++)
	{
		for (int i = 0; i <= (k-1); i++)
		{
			sum = 0;
			for (int j = 0; j <= (i-1); j++) sum = sum + a[i][j] * a[k][j];
			a[k][i] = (a[k][i] - sum) / a[i][i];
		}
		sum = 0;
		for (int j = 0; j <= (k-1); j++) sum = sum + pow(a[k][j],2);
		a[k][k] = sqrt(a[k][k] - sum);
	}
	copymatrix(a,b,n,n);
}

int main()
{
	double A[n][n] = {{6,15,55},{15,55,225},{55,225,979}}; // Matrix A, Initialize
	cout << "A = \n"; printmatrix(A,n,n,0); // Matrix A, Print
	double mA[n][n]; // Matrix mA, Intialize
	cLudecomp(A,mA); // mA ~ A
	cout << "mA = \n"; printmatrix(mA,n,n,1); // Matrix mA, Print
	
	cout << "\n";
	system("pause");
	return 0;
}