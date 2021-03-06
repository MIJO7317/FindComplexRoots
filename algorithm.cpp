#include "fparser.hh"
#include<iostream>
#include<algorithm>
#include<complex>
#include<cmath>
#include<vector>

#define pi long double(4*atan(1))
#define e long double(exp(1))

FunctionParser_cld fparser;
std::string function;

namespace FindCR
{
	#define pi long double(4*atan(1))

	typedef std::complex<long double> complex;

	complex f(complex z)
	{
		return fparser.Eval(&z);
	}

	long double AngleBetween(complex z1, complex z2)
	{
		long double angle = fabsl(arg(z1) - arg(z2));
		return angle > pi ? 2 * pi - angle : angle;
	}

	int sign(long double x)
	{
		return x == 0 ? 0 : x<0 ? -1 : 1;
	}

	void FindRoots(long double radius_square, long double epsilon, long double x0, long double y0, std::vector<complex> &roots)
	{
		complex z_prev(0,0), z_curr;
		long double sum = 0, d_alpha = 0.01;;
		for (long double alpha = 0; alpha < 2 * pi; alpha += d_alpha)
		{
			long double real = x0 + radius_square * cosl(alpha) / std::fmaxl(fabsl(sinl(alpha)), fabsl(cosl(alpha)));
			long double imag = y0 + radius_square * sinl(alpha) / std::fmaxl(fabsl(sinl(alpha)), fabsl(cosl(alpha)));
			z_curr = complex(real,imag);
			z_curr = f(z_curr);
			long double angle = AngleBetween(z_curr, z_prev);
			int sign_sum = sign(std::real(z_curr) * std::imag(z_prev) - std::imag(z_curr) * std::real(z_prev));
			sum += sign_sum*angle;
			z_prev = complex(std::real(z_curr), std::imag(z_curr));
		}

		if (fabsl(sum) / 2 / pi >= 1-d_alpha)
		{
			if (radius_square < epsilon)
			{
				bool is_find = false;
				for (auto el : roots)
				{
					if (fabsl(real(el) - x0) <= 2*epsilon && fabsl(imag(el) - y0) <= 2*epsilon)
					{
						is_find = true;
						break;
					}
				}
				if (!is_find)
					roots.push_back(complex(round(x0/epsilon)*epsilon,round(y0/epsilon)*epsilon));
			}
			else
			{
				FindRoots(radius_square / 2 + radius_square / 1000, epsilon, x0 + radius_square / 2, y0 + radius_square / 2, roots);
				FindRoots(radius_square / 2 + radius_square / 1000, epsilon, x0 + radius_square / 2, y0 - radius_square / 2, roots);
				FindRoots(radius_square / 2 + radius_square / 1000, epsilon, x0 - radius_square / 2, y0 + radius_square / 2, roots);
				FindRoots(radius_square / 2 + radius_square / 1000, epsilon, x0 - radius_square / 2, y0 - radius_square / 2, roots);
			}
		}
	}

}

std::complex<long double> a, b;
std::vector<std::complex<long double>> roots;
long double epsilon, radius_square, x_centre, y_centre;

int main()
{
	fparser.AddConstant("pi", pi);
	fparser.AddConstant("e", e);
	while (true)
	{
		std::cout << "f(z) = ";
		std::getline(std::cin, function);
		if (std::cin.fail()) return 0;

		int res = fparser.Parse(function, "z");
		if (res < 0) break;

		std::cout << std::string(res + 7, ' ') << "^\n"
			<< fparser.ErrorMsg() << "\n\n";
	}
	a = std::complex<long double>(1, 0);
	b = std::complex<long double>(1, 1);
	std::cout << "Epsilon = ";
	std::cin >> epsilon;
	std::cout << "Radius square = ";
	std::cin >> radius_square;
	std::cout << "Centre square = ";
	std::cin >> x_centre >> y_centre;
	FindCR::FindRoots(radius_square, epsilon, x_centre, y_centre, roots);
	std::cout << "Roots:\n";
	for (auto el : roots)
	{
		std::cout.precision(std::to_string(int(1/epsilon)).length());
		std::cout << el << std::endl;
	}
}
