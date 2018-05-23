// 03/05/2014 - ENGR 2450 - Meine, Joel
// Problems 17.20, 20.22

// Polynomial Regression

#include <iostream>
#include <iomanip>
#include <math.h>
using namespace std;

const int nn = 20;
const int N = 7; // Number of Data Points
const int M2 = 2; // Polynomial, Second Order
const int M3 = 3; // Polynomial, Third Order
const int M2p1 = M2 + 1;  // Polynomial, Second Order plus One
const int M3p1 = M3 + 1;  // Polynomial, Third Order plus One

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

void NLRegress(double Z[][nn],double Y[],double A[],int n,int m,int p,double tol,int er)
{
	if (p == 2) // Polynomial, Second Order
	{
		double ZT[M2p1][nn]; double ZTZ[M2p1][nn]; double ZTZI[M2p1][nn]; double ZTY[M2p1];
		mtranspose(Z,ZT,n,(m+1));
		multiply_matrices(ZT,Z,ZTZ,(m+1),(m+1),n);
		MatrixInverse(ZTZ,ZTZI,(m+1),tol,er);
		multiply_matrix_to_vector(ZT,Y,ZTY,(m+1),n);
		multiply_matrix_to_vector(ZTZI,ZTY,A,(m+1),(m+1));
	}
	else if (p == 3) // Polynomial, Third Order
	{
		double ZT[M3p1][nn]; double ZTZ[M3p1][nn]; double ZTZI[M3p1][nn]; double ZTY[M3p1];
		mtranspose(Z,ZT,n,(m+1));
		multiply_matrices(ZT,Z,ZTZ,(m+1),(m+1),n);
		MatrixInverse(ZTZ,ZTZI,(m+1),tol,er);
		multiply_matrix_to_vector(ZT,Y,ZTY,(m+1),n);
		multiply_matrix_to_vector(ZTZI,ZTY,A,(m+1),(m+1));
	}
}

void BuildZP(double x[],int n,int m,double Z[][nn])
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < (m+1); j++)
			Z[i][j] = pow(x[i],j);
	}
}

int main()
{
	// Problem 17.20 - Polynomial, Second Order
	double x[N] = {0.2,0.5,0.8,1.2,1.7,2,2.3};
	double y[N] = {500,700,1000,1200,2200,2650,3750};

	double Z2[N][nn];
	BuildZP(x,N,M2,Z2);

	double a2[M2p1];
	NLRegress(Z2,y,a2,N,M2,2,0.0000000001,0);

	double yf[N];
	multiply_matrix_to_vector(Z2,a2,yf,N,M2p1);

	std::cout << "Chapter 17 - Problem 17.20" << std::endl;
	std::cout << "============================================" << std::endl;
	std::cout << "Data Fitting: Polynomial, Second Order" << std::endl;
	cout << "\n";
	std::cout << "yf = a0 + a1*x + a2*x^2" << std::endl;
	cout << "\n";
	std::cout << "a0 = " << setprecision(6) << a2[0] << std::endl;
	std::cout << "a1 = " << setprecision(6) << a2[1] << std::endl;
	std::cout << "a2 = " << setprecision(6) << a2[2] << std::endl;
	std::cout << "********************************************" << std::endl;
	std::cout << " x   y      yf" << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	for (int i = 0; i < N; i++)
	{
		cout << setw(3) << x[i];
		cout << setw(6) << y[i];
		cout << setw(9) << setprecision(5) << yf[i] << endl;
	}
	std::cout << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	cout << "\n";
	
	// Problem 20.22 - Polynomial, Third Order
	double T[N] = {0,5,10,15,20,25,30}; // Temperature (C)
	double DO[N] = {12.9,11.3,10.1,9.03,8.17,7.46,6.85}; // Dissolved Oxygen (mg/L)

	double Z3[N][nn];
	BuildZP(T,N,M3p1,Z3);

	double a3[M3p1];
	NLRegress(Z3,DO,a3,N,M3,3,0.0000000001,0);

	double DOf[N];
	multiply_matrix_to_vector(Z3,a3,DOf,N,M3p1);

	int t = 8; // Temperature = 8 C
	double DOe = a3[0] +a3[1]*t + a3[2]*pow(t,2) + a3[3]*pow(t,3); // DO(T)

	std::cout << "Chapter 20 - Problem 20.22" << std::endl;
	std::cout << "============================================" << std::endl;
	std::cout << "Data Fitting: Polynomial, Third Order" << std::endl;
	cout << "\n";
	std::cout << "DOf = a0 + a1*T + a2*T^2 + a3*T^3" << std::endl;
	cout << "\n";
	std::cout << "a0 = " << setprecision(6) << a3[0] << std::endl;
	std::cout << "a1 = " << setprecision(6) << a3[1] << std::endl;
	std::cout << "a2 = " << setprecision(6) << a3[2] << std::endl;
	std::cout << "a3 = " << setprecision(6) << a3[3] << std::endl;
	cout << "\n";
	std::cout << "Dissolved Oxygen, Fitted (mg/L), DOf" << std::endl;
	std::cout << "Dissolved Oxygen, Recorded (mg/L), DO" << std::endl;
	std::cout << "Temperature (C), T" << std::endl;
	std::cout << "Concentration of Chloride (g/L), c = 10" << std::endl;
	cout << "\n";
	std::cout << "DOf(T=" << t << ") = " << DOe << std::endl;
	std::cout << "********************************************" << std::endl;
	std::cout << " T   DO     DOf" << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	for (int i = 0; i < N; i++)
	{
		cout << setw(3) << T[i];
		cout << setw(6) << DO[i];
		cout << setw(9) << setprecision(5) << DOf[i] << endl;
	}
	std::cout << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	
	cout << "\n";
	system("pause");
	return 0;
}