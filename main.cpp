#include <iostream>
#include <stdarg.h>
#include <math.h>

class Polynomial
{
	// p_0 * x^n + p_1 * x^(n-1) + ... p_(n-1)*x + p_n
public:
	long double* coefficients_list;
	int degree;

	Polynomial(int n, ...)
	{
		degree = n;
		coefficients_list = new long double[n + 1];
		// указатель на начало обрабатываемых аргументов
		va_list coefficients;
		// начало обработки аргументов
		va_start(coefficients, n);
		for (int i = 0; i <= n; i++)
		{
			coefficients_list[i] = va_arg(coefficients, long double);
		}
		va_end(coefficients);
	}

	explicit Polynomial(int n, long double* coeficients)
	{
		coefficients_list = coeficients;
		degree = n;
	}

	~Polynomial()
	{
		delete[] coefficients_list;
	}

	void printCoefs()
	{
		for (int i = 0; i <= degree; i++)
		{
			std::cout << coefficients_list[i] << std::endl;
		}
	}

	Polynomial* calcDerivative() const
	{
		long double* derivative_coeficients = new long double[degree];
		for (int i = 0; i < degree; i++)
		{
			long double n_i = degree - i;
			long double p_i = coefficients_list[i];
			derivative_coeficients[i] = p_i * n_i;
		}
		Polynomial* derivative = new Polynomial(degree - 1, derivative_coeficients);
		return derivative;
	}

	long double calcF(long double x)
	{
		long double result = 0;
		for (int i = 0; i <= degree; i++)
		{
			long double x_i = pow(x, i);
			long double p_i = coefficients_list[degree - i];
			result += x_i * p_i;
		}
		return result;
	}
};


long double calcAverageX(long double x_1, long double x_2)
{
	return (x_1 + x_2) / 2;
}

long double findXL(Polynomial* f, long double x_R)
{
	long double delta = 1;
	long double x_L = x_R - delta;
	long double f_x = f->calcF(x_R - delta);
	Polynomial* df = f->calcDerivative();
	long double df_x = df->calcF(x_R - delta);
	while ((f_x > 0  && df_x > 0) || (f_x < 0 && df_x < 0))
	{
		x_L -= delta;
		f_x = f->calcF(x_L);
		df_x = df->calcF(x_L);
	}
	return x_L;
}

// x_L < x_R
long double bisection(Polynomial* f, long double x_L, long double x_R, long double epsilon)
{
	// Ищем левую границу отрезка для бисекции
	if (x_L == INFINITY)
	{
		x_L = findXL(f, x_R);
	}

	// Ищем правую границу для бисекции
	else if (x_R == INFINITY)
	{
		long double delta = 10;
		while (abs(f->calcF(x_L + delta)) >= epsilon)
		{
			delta *= 2;
			x_R = delta;
		}
	}

	long double x_average = calcAverageX(x_L, x_R);
	long double f_x = f->calcF(x_average);
	while (abs(f_x) >= epsilon)
	{
		// Ищем решение на [x_average; x_R]
		if (f_x < -epsilon)
		{
			x_L = x_average;
			x_average = calcAverageX(x_average, x_R);
		}
			// Ищем решение на [x_L; x_average]
		else
		{
			x_R = x_average;
			x_average = calcAverageX(x_L, x_average);
		}
		f_x = f->calcF(x_average);
	}
	return x_average;
}

long double calcD(long double a, long double b, long double c)
{
	return pow(b, 2) - 4 * a * c;
}

long double* findRoots(Polynomial* f, long double epsilon)
{
	long double* roots;
	// Полином первой степени
	if (f->degree == 1)
	{
		roots = new long double;
		if (f->coefficients_list[0] == 0)
		{
			return nullptr;
		}
		else if (((f->calcF(0) < -epsilon) && (f->coefficients_list[0] < 0)) ||
				 ((f->calcF(0) > epsilon) && (f->coefficients_list[0] > 0)))
		{
			roots[0] = bisection(f, INFINITY, 0.0, epsilon);
			return roots;
		}
		else
		{
			roots[0] = bisection(f, 0.0, INFINITY, epsilon);
			return roots;
		}
	}

	// Полином второй степени
	else if (f->degree == 2)
	{
		long double D = calcD(f->coefficients_list[0], f->coefficients_list[1], f->coefficients_list[2]);

		// Определяем количество корней

		// Один корень
		if (D == 0)
		{
			roots = new long double;
			Polynomial *d_f = f->calcDerivative();
			roots[0] = *findRoots(d_f, epsilon);
			return roots;
		}
		else if (D > 0)
		{

		}
		else
		{

		}
	}

	// Полином третьей степени
	else if (f->degree == 3)
	{
		Polynomial* derivative = f->calcDerivative();
		long double c = derivative->coefficients_list[2];

		// Определяем количество корней

		long double D = calcD(derivative->coefficients_list[0], derivative->coefficients_list[1],
				derivative->coefficients_list[2]);
		// Один корень. f'(x) = 0 или f'(x) > 0
		if ((D < -epsilon) || (abs(D) <= epsilon))
		{
			roots = new long double;
			if (abs(c) < epsilon)
			{
				*roots = 0;
				return roots;
			}
			else
			{
				if (c > epsilon)
				{
					*roots = bisection(derivative, c, INFINITY, epsilon);
					return roots;
				}
				else
				{
					*roots = bisection(derivative, INFINITY, c, epsilon);
					return roots;
				}
			}
		}

			// Два корня.
		else
		{

		}
	}

}

int main()
{
	int n = 0;
	std::cout << "Введите степень полинома" << std::endl;
	std::cin >> n;
	auto* coeficients = new long double[n + 1];
	std::cout << "Введите коэффициенты многочлена, начиная с коэффициента при наибольшей степени" << std::endl;
	for (int i = 0; i <= n; i++)
	{
		std::cin >> coeficients[i];
	}

	Polynomial* f = new Polynomial(n, coeficients);
	std::cout << std::endl;
	//long double root_x = bisection(f, -10, 10, 0.05);
	long double root_x = *findRoots(f, 0.1);
	std::cout << "root = " << root_x << std::endl;
	Polynomial* derivative = f->calcDerivative();
	derivative->printCoefs();
	delete[] coeficients;
	return 0;
}
