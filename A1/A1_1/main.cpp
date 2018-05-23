// 01/21/2014 - ENGR 2450 - Meine, Joel
// Chapter 1 - Problem 1.8

#include <iostream>
#include <math.h>
using namespace std;

int main()
{
	double A = 1200; // Area (m^2)
	double Q = 500; // Flow (m^3/d)
	double a = 300; // Constant
	double y = 0; // Depth_initial (m)
	
	std::cout << "Chapter 1 - Problem 1.8" << std::endl;
	std::cout << "============================" << std::endl;
	std::cout << "Area, A(m^2) = " << A << std::endl;
	std::cout << "Flow, Q(m^3/d) = " << Q << std::endl;
	std::cout << "Constant, a = " << a << std::endl;
	std::cout << "Depth_initial, y(m) = " << y << std::endl;
	std::cout << "----------------------------" << std::endl;
	std::cout << " t(s)     y(m)              " << std::endl;
	std::cout << "----------------------------" << std::endl;

	double  t_f = 10,  t_s = 0.5;
	int t = 0, T = t_f / t_s;

	printf(" %2.1f %12.10f \n", t, y);
	for (t; t < T; t++)
	{
		// value_new = value_old + step_size*slope
		y = y + t_s*(3*(Q/A)*pow(sin(t*t_s),2) - (a*pow((1+y),1.5))/A);
		printf(" %2.1f %12.10f \n", t*t_s + t_s, y);
	}
	std::cout << "----------------------------" << std::endl;
	system("pause");
	return 0;
}