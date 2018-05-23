// 03/25/2014 - ENGR 2450 - Meine, Joel
// Problems 21.9, 24.22

// Trapezoidal Rule (Multiple-Segment)

#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
using namespace std;

double Trapm(double h,int n,vector<double> f)
{
	double sum = f[0];
	for (int i = 1; i < n; i++)
		sum = sum + 2 * f[i];
	sum = sum + f[n];
	return(h * (sum/2));
}

double V(double t)
{
	const double g = 9.8; // Gravitational Constant (m/s^2), g
	const double m = 68.1; // Mass of Object (kg), m
	const double cd = 0.25; // Drag Coefficient (kg/m), cd
	double v = sqrt((g*m)/cd)*tanh(sqrt((g*cd)/m)*t); // Velocity of Object (m/s), v
	return(v);
}

void funcEval(double a,double b,int n,double& h,double& I)
{
	h = (b-a)/n;
	vector<double> y;
	for (int i = 0; i <= n; i++)
	{
		y.push_back(V(a+i*h));
	}
	I = Trapm(h,n,y);
	return;
}

void dataEval(double a,double b,int n,vector<double> y,double& I)
{
	double h = (b-a)/n;
	I = Trapm(h,n,y);
	return;
}

int main()
{
	double h = 0; // Size of Subinterval, h
	
	// Problem 21.9 - Integral of a Function
	const double ti = 0; // Initial Time (s), t
	const double tf = 10; // Final Time (s), t
	double s = 0; // Position of Object (m), s

	std::cout << "Chapter 21 - Problem 21.9" << std::endl;
	std::cout << "============================================" << std::endl;
	std::cout << "Initial Time (s), ti = " << setprecision(0) << ti << std::endl;
	std::cout << "Final Time (s), tf = " << setprecision(0) << tf << std::endl;
	std::cout << "Number of Subintervals, n" << std::endl;
	std::cout << "Size of Subinterval, h" << std::endl;
	std::cout << "Position of Object (m), s" << std::endl;
	std::cout << "********************************************" << std::endl;
	std::cout << " n     h     s" << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	const int n_size = 6;
	int n[] = {5,10,50,100,500,1000};
	for (int i = 0; i < n_size; i++)
	{
		funcEval(ti,tf,n[i],h,s);
		cout << setw(5) << n[i];
		cout << setw(6) << setprecision(3) << h;
		cout << setw(13) << setprecision(10) << s << endl;
	}
	std::cout << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	cout << "\n";

	// Problem 24.22 - Integral for Tabulated Data
	vector<double> l = {0,30,60,90,120,150,180,210,240}; // Height (m), l
	vector<double> F = {0,340,1200,1600,2700,3100,3200,3500,3800}; // Force (N/m), F(l)
	vector<double> lF;
	for (int i = 0; i < l.size(); i++)
		lF.push_back(l[i]*F[i]);
	double li = l.front(); // Initial Height (m), li
	double lf = l.back(); // Final Height (m), lf
	double T = 0; // Total Force (N), T
	double d = 0; // Line of Action (m), d

	std::cout << "Chapter 24 - Problem 24.22" << std::endl;
	std::cout << "============================================" << std::endl;
	std::cout << "Initial Height (m), li = " << setprecision(0) << li << std::endl;
	std::cout << "Final Height (m), lf = " << setprecision(0) << lf << std::endl;
	std::cout << "********************************************" << std::endl;
	
	dataEval(li,lf,F.size()-1,F,T);
	std::cout << "Total Force (N), T = " << setprecision(10) << T << std::endl;
	
	double dT = 0;
	dataEval(li,lf,lF.size()-1,lF,dT);
	d = dT/T;
	std::cout << "Line of Action (m), d = " << setprecision(10) << d << std::endl;
	std::cout << "++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	cout << "\n";

	system("pause");
	return 0;
}