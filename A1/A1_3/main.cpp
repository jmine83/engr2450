// 01/21/2014 - ENGR 2450 - Meine, Joel
// Chapter 3 - Problem 3.13

#include <iostream>
#include <math.h>
using namespace std;

int main()
{
	double x_guess = 5.0; // Guess_initial
	int n = 5; // Decimals_number
	double es = 0.5*pow(10,2-n); // Error Criteria
	int maxit = 20; // Maximum Number of Iterations
	double sol; // Solution_current
	double solold; // Solution_previous
	int iter; // Iteration
	double a = 122.5; // Square Root Input
	double ea; // Relative Error

	std::cout << "Chapter 3 - Problem 3.13" << std::endl;
	std::cout << "============================" << std::endl;
	std::cout << "Guess_initial, x_guess = " << x_guess << std::endl;
	std::cout << "Decimals_number, n = " << n << std::endl;
	std::cout << "Error Criteria, es = " << es << std::endl;
	std::cout << "Maximum Number of Iterations, maxit = " << maxit << std::endl;
	std::cout << "Square Root Input, a = " << a << std::endl;
	std::cout << "----------------------------" << std::endl;
	std::cout << "iter.  x_current  rel. error" << std::endl;
	std::cout << "----------------------------" << std::endl;

	iter = 1; sol = x_guess; ea = 100;
	printf(" %2i %12.10f %12.10f ", iter, sol, ea);
	do {
		solold = sol;
		sol = 0.5*(sol+(a/sol)); // Divide-and-Average Method
		iter = iter + 1;
		if (sol != 0) ea = abs((sol-solold)/sol)*100;
		printf("\n %2i %12.10f %12.10f ", iter, sol, ea);
	} while (ea > es && iter < maxit);
	if (ea <= es)
	{
		std::cout << "\n----------------------------" << std::endl;
		printf("x_solution = %12.10f \n", sol);
		printf("\nrel. error = %12.10f \n", ea);
		std::cout << "----------------------------" << std::endl;
	}
	else if (iter >= maxit)
		std::cout << "No Solution due to Divergence" << std::endl;
	system("pause");
	return 0;
}