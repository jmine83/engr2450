// 02/05/2014 - ENGR 2450 - Meine, Joel
// Chapter 6 - Problem 6.30

#include <iostream>
#include <math.h>
using namespace std;

// Spherical Tank Volume
const double R = 3;
const double PI = 3.141592653589793;

double V(double h,double v)
{
	double Vr = PI * pow(h,2) * ((3*R-h)/3) - v;
	return Vr;
}

double dV(double dh)
{
	double dv = PI*2*R*dh - PI*pow(dh,2);
	return dv;
}

// Problem Function
double mainF(double xm,double ym)
{
	double Y = V(xm,ym);
	return Y;
}

double maindF(double dxm)
{
	double dY = dV(dxm);
	return dY;
}

void mainR(double i0,double i1,double i2,double i3,double i4,int i5)
{
	printf("%2.0f %2.0f  %2.3f %2.4f  %2.4f  %2i \n",i0,i1,i2,i3,i4,i5);
}

// Test Function
double testF(double xt)
{
	double G = pow(xt,2) - 5*xt + 6;
	return G;
}

double testdF(double dxt)
{
	double dG = 2*dxt - 5;
	return dG;
}

void testR(double j0,double j1,double j2,double j3,int j4)
{
	printf("%2.1f  %2.4f %2.4f %2.4f %2i \n",j0,j1,j2,j3,j4);
}

// Newton-Raphson Method
void NewtonRaphson(int F,double C,double X0,double Es,int Imax,int Ilim)
{
	int iter = 0; // Iteration
	double ea = Es + 1; // Absolute Error, actual (initial)
	double xr = X0; // Root of Function (initial)
	double xrold = 0; // Root of Function (previous)
	double fr = 0; // Function Value at Root (initial)

	do
	{
		xrold = xr;
		if (F == 0) xr = xrold - (mainF(xrold,C)/maindF(xrold));
		else if (F == 1) xr = xrold - (testF(xrold)/testdF(xrold));
		iter++;
		if (xr != 0) ea = abs((xr-xrold)/xr)*100;
		if (F == 0 && iter <= Ilim)
		{
			fr = mainF(xr,C);
			printf(" %2i    %2.1f  %2.2f  %2.3f  %2.2f  %2.2f \n",iter,C,X0,xr,fr,ea);
		}
	} while (ea >= Es && iter < Imax);
	if (ea < Es)
	{
		if (F == 0 && Ilim <= -1)
		{
			fr = mainF(xr,C);
			mainR(C,X0,xr,fr,ea,iter);
		}
		else if (F == 1)
		{
			fr = testF(xr);
			testR(X0,xr,fr,ea,iter);
		}
	}
	else if (iter >= Imax) std::cout << "No Solution due to Divergence" << std::endl;
}

int main()
{
	const double es = 0.0005; // Error Criteria;
	const int imax = 100; // Maximum Number of Iterations

	std::cout << "Chapter 6 - Problem 6.30" << std::endl;
	std::cout << "============================================" << std::endl;
	std::cout << "Guess_initial, x0" << std::endl;
	std::cout << "Root of Function, xr" << std::endl;
	std::cout << "Function Value at Root, f(xr)" << std::endl;
	std::cout << "Error Criteria, es = " << es << std::endl;
	std::cout << "Actual Error, ea" << std::endl;
	std::cout << "Maximum Number of Iterations, imax = " << imax << std::endl;
	std::cout << "Iterations, iter" << std::endl;
	std::cout << "********************************************" << std::endl;

	std::cout << "Test Function" << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	std::cout << "x0   xr     f(xr)  ea     iter" << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	NewtonRaphson(1,0,5.0,es,imax,0);
	NewtonRaphson(1,0,1.0,es,imax,0);
	std::cout << "********************************************" << std::endl;

	std::cout << "Problem Function" << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	std::cout << "Volume of Spherical Tank, V(m^3)" << std::endl;
	std::cout << "Spherical Tank Radius, R(m) = " << R << std::endl;
	std::cout << "Depth of Liquid in Tank, h(m)" << std::endl;
	std::cout << "Depth_initial, h0(m)" << std::endl;
	std::cout << "Depth_root, hr(m)" << std::endl;
	std::cout << "Function Value at Depth_root, f(hr)" << std::endl;
	std::cout << "Iteration, iter." << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	std::cout << "iter.  V     h0     hr      f(hr)  ea" << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	NewtonRaphson(0,30,10,es,imax,3);
	std::cout << "********************************************" << std::endl;
	std::cout << "V  h0  hr     f(hr)   ea      iter" << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	NewtonRaphson(0,30,10,es,imax,-1);
	std::cout << "******************************************** \n" << std::endl;

	system("pause");
	return 0;
}