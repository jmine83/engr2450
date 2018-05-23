// 03/05/2014 - ENGR 2450 - Meine, Joel
// Problem 17.18

// Multiple-Linear Regression

#include <iostream>
#include <iomanip>
#include <math.h>
using namespace std;

const int nn = 20;
const int N = 9; // Number of Data Points
const int M = 2; // Independent Variable Set, M Order
const int Mp1 = M + 1; // Independent Variable Set, M plus One

void mtranspose(double M[][nn],double MT[][nn],int n,int m)
{//Calculates transpose of matrix M(nxm) as MT(mxn)
	int i; int j;
	for (i = 0; i < n; i++){
		for (j = 0; j < m; j++){
			MT[j][i] = M[i][j];
		};
	};
}

void multiply_matrices(double A[][nn],double B[][nn],double C[][nn],\
	int n,int m,int p)
{//calculates C(nxm) = A(nxp)*B(pxm)
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

void multiply_matrix_to_vector(double A[][nn],double u[],double v[],\
	int n,int m)
{//calculates v(n) = A(nxm)*u(m)
	int i; int k;
	for (i = 0; i < n; i++){
		v[i] = 0;
		for (k = 0; k < m; k++){
			v[i] = v[i] + A[i][k] * u[k];
		};
	};
}

void Pivot(double A[][nn],int O[],double s[],int n,int k)
{//Pivot routine for LU decomposition 
	//Use it together with the MatrixInverse subroutine
	int p; double big; double dummyB; int ii; int dummyP;
	p = k; big = abs(A[O[k]][k] / s[O[k]]);
	for (ii = k + 1; ii<n; ii++)
	{
		dummyB = abs(A[O[ii]][k] / s[O[ii]]);
		if (dummyB>big)
		{
			big = dummyB; p = ii;
		};
	};
	dummyP = O[p]; O[p] = O[k]; O[k] = dummyP;
};

void Substitute(double A[][nn],int O[],int n,double b[],double x[])
{//Substitution routine for LU decomposition
	//Use it together with the MatrixInverse subroutine
	double sum; int i; int j;
	for (i = 1; i<n; i++)
	{
		sum = b[O[i]];
		for (j = 0; j <= i - 1; j++)
		{
			sum = sum - A[O[i]][j] * b[O[j]];
		};
		b[O[i]] = sum;
	};
	x[n - 1] = b[O[n - 1]] / A[O[n - 1]][n - 1];
	for (i = n - 2; i >= 0; i--)
	{
		sum = 0.0;
		for (j = i + 1; j<n; j++)
		{
			sum = sum + A[O[i]][j] * x[j];
		};
		x[i] = (b[O[i]] - sum) / A[O[i]][i];
	};
};

void Decompose(double A[][nn],int n,double tol,int O[],double s[],int er)
{//Decomposition step in the LU decomposition algorithm
	//Use it together with the MatrixInverse subroutine
	int i; int j; int k; double factor;
	for (i = 0; i<n; i++)
	{
		O[i] = i; s[i] = abs(A[i][0]);
		for (j = 1; j<n; j++)
		{
			if (abs(A[i][j])>s[i])
			{
				s[i] = abs(A[i][j]);
			};
		};
	};
	for (k = 0; k<n - 1; k++)
	{
		Pivot(A,O,s,n,k);
		if (abs(A[O[k]][k] / s[O[k]])<tol)
		{
			er = -1; cout << "trouble" << endl;
			break;
		};
		for (i = k + 1; i<n; i++)
		{
			factor = A[O[i]][k] / A[O[k]][k];
			A[O[i]][k] = factor;
			for (j = k + 1; j<n; j++)
			{
				A[O[i]][j] = A[O[i]][j] - factor*A[O[k]][j];
			}
		};
	};
	if (abs(A[O[k]][k] / s[O[k]]) < tol)
	{
		er = -1; cout << "trouble" << endl;
	}
};

void MatrixInverse(double A[][nn],double AI[][nn],int n,double tol,int er)
{//Calculates the inverse of matrix A(nxn), i.e., AI(nxn)
	//This subroutine is the main driver for Matrix Inverse with LU decomposition algorithm
	int O[nn]; double s[nn]; double b[nn]; double x[nn];
	int i; int j;
	Decompose(A,n,tol,O,s,er);
	if (er == 0)
	{
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				if (i == j)
				{
					b[j] = 1;
				}
				else
				{
					b[j] = 0;
				}
			}
			Substitute(A,O,n,b,x);
			for (j = 0; j < n; j++)
			{
				AI[j][i] = x[j];
			}
		}
	}
	else
	{
		cout << "Matrix is singular" << endl;
	}
}

void NLRegress(double Z[][nn],double Y[],double A[],int n,int m,double tol,int er)
{
	double ZT[Mp1][nn]; double ZTZ[Mp1][nn]; double ZTZI[Mp1][nn]; double ZTY[Mp1];
	mtranspose(Z,ZT,n,(m+1));
	multiply_matrices(ZT,Z,ZTZ,(m+1),(m+1),n);
	MatrixInverse(ZTZ,ZTZI,(m+1),tol,er);
	multiply_matrix_to_vector(ZT,Y,ZTY,(m+1),n);
	multiply_matrix_to_vector(ZTZI,ZTY,A,(m+1),(m+1));
}

void BuildZM(double X[][nn],int n,int m,double Z[][nn])
{
	for (int i = 0; i < n; i++)
	{
		Z[i][0] = 1.0;
		for (int j = 1; j < (m+1); j++)
			Z[i][j] = X[i][j-1];
	}
}

int main()
{
	// Problem 17.18 - Multiple Linear Regression
	double x[N][nn] = {{0,0},{1,1},{1,2},{2,1},{2,2},{3,1},{3,2},{4,1},{4,2}};
	double y[N] = {15.1,17.9,12.7,25.6,20.5,35.1,29.7,45.4,40.2};

	double Z[N][nn];
	BuildZM(x,N,M,Z);

	double a[Mp1];
	NLRegress(Z,y,a,N,M,0.0001,0);

	double yf[N];
	multiply_matrix_to_vector(Z,a,yf,N,Mp1);

	std::cout << "Chapter 17 - Problem 17.18" << std::endl;
	std::cout << "============================================" << std::endl;
	std::cout << "yf = a0 + a1*x1 + a2*x2" << std::endl;
	cout << "\n";
	std::cout << "Slope for x2, a2 = " << setprecision(4) << a[2] << std::endl;
	std::cout << "Slope for x1, a1 = " << setprecision(4) << a[1] << std::endl;
	std::cout << "Intercept, a0 = " << setprecision(4) << a[0] << std::endl;
	std::cout << "********************************************" << std::endl;
	std::cout << " x1  x2   y      yf" << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	for (int i = 0; i < N; i++)
	{
		cout << setw(3) << x[i][0];
		cout << setw(4) << x[i][1];
		cout << setw(7) << y[i];
		cout << setw(9) << setprecision(5) << yf[i] << endl;
	}
	std::cout << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	
	cout << "\n";
	system("pause");
	return 0;
}