// 01/21/2014 - ENGR 2450 - Meine, Joel
// Chapter 2 - Problem 2.22

#include <iostream>
#include <math.h>
using namespace std;

int main()
{
	double x_i = 0; // Position_initial (ft)
	double x_f = 10; // Position_final (ft)
	double dx = 0.5; // Position_step (ft)
	double x = 0; // Position_current (ft)
	
	std::cout << "Chapter 2 - Problem 2.22" << std::endl;
	std::cout << "============================" << std::endl;
	std::cout << "Position_initial, x_i(ft) = " << x_i << std::endl;
	std::cout << "Position_final, x_f(ft) = " << x_f << std::endl;
	std::cout << "Position_step, dx(ft) = " << dx << std::endl;
	std::cout << "Position_current, x(ft)" << std::endl;
	std::cout << "Beam Displacement, u_y(ft)" << std::endl;
	std::cout << "----------------------------" << std::endl;
	std::cout << "x(ft)   u_y(ft)             " << std::endl;
	std::cout << "----------------------------" << std::endl;

	int n = (x_f-x_i)/dx + 1;
	for (int i = 1; i <= n; i++)
	{
		x = x_i + (i - 1)*dx;
		double u_y = 0;
		if (x <= 0)
			u_y = (57.0/6.0)*pow(x,3) - 238.25*x;
		else if (x > 0 && x <= 5.0)
			u_y = (-5.0/6.0)*pow(x, 4) + (57.0/6.0)*pow(x,3) - 238.25*x;
		else if (x > 5.0 && x <= 7.0)
			u_y = (-5.0/6.0)*(pow(x,4)-pow(x-5,4)) + (57.0/6.0)*pow(x,3) - 238.25*x;
		else if (x > 7.0 && x <= 8.0)
			u_y = (-5.0/6.0)*(pow(x,4)-pow(x-5,4)) + 75.0*pow(x-7,2) + (57.0/6.0)*pow(x,3) - 238.25*x;
		else if (x > 8.0)
			u_y = (-5.0/6.0)*(pow(x,4) - pow(x-5,4)) + (15.0/6.0)*pow(x-8,3) + 75.0*pow(x-7,2) + (57.0/6.0)*pow(x,3) - 238.25*x;
		else
			std::cout << "ERROR!" << std::endl;
		printf(" %2.1f %12.10f \n", x, u_y);
	}
	std::cout << "----------------------------" << std::endl;
	system("pause");
	return 0;
}