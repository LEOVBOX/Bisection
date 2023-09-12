#include <iostream>
#include <stdarg.h>
#include <math.h>

class Polynomial
{
	// p_0 * x^n + p_1 * x^(n-1) + ... p_(n-1)*x + p_n
	// p_0 > 0
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
	long double f_x = f->calcF(x_L);
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

long double findXR(Polynomial* f, long double x_L)
{
	long double delta = 1;
	long double x_R = x_L + delta;
	long double f_x = f->calcF(x_R);
	Polynomial* df = f->calcDerivative();
	long double df_x = df->calcF(x_R - delta);
	while ((f_x > 0 && df_x < 0) || (f_x < 0 && df_x > 0))
	{
		x_R += delta;
		f_x = f->calcF(x_R);
		df_x = df->calcF(x_R);
	}
	return x_R;
}

long double bisection(Polynomial* f, long double x_L, long double x_R, long double epsilon)
{
	// x_L < x_R
	// Ищем левую границу отрезка для бисекции
	if (x_L == INFINITY)
	{
		x_L = findXL(f, x_R);
	}

	// Ищем правую границу для бисекции
	else if (x_R == INFINITY)
	{
		x_R = findXR(f, x_L);
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
		roots = new long double[2];
		long double k = f->coefficients_list[0];
		long double b = f->coefficients_list[1];
		if (k == 0)
		{
			return nullptr;
		}
		else
		{
			roots[0] = 1;
			roots[1] = -(b/k);
		}
	}

	// Полином второй степени
	else if (f->degree == 2)
	{
		long double a = f->coefficients_list[0];
		long double b = f->coefficients_list[1];
		long double c = f->coefficients_list[2];

		long double D = calcD(a, b, c);

		// Определяем количество корней

		// Один корень
		if (D == 0)
		{
			roots = new long double[2];
			roots[0] = 1;
			roots[1] = -(b/(2*a));
		}
		else if (D > 0)
		{
			roots = new long double[3];
			roots[0] = 2;
			roots[1] = (-b - sqrt(D))/(2*a);
			roots[2] = (-b + sqrt(D))/(2*a);
			return roots;
		}
		else
		{
			return nullptr;
		}
	}

	// Полином третьей степени
	else
	{
		Polynomial* derivative = f->calcDerivative();
		long double c = f->coefficients_list[3];

		long double D_df = calcD(derivative->coefficients_list[0], derivative->coefficients_list[1],
				derivative->coefficients_list[2]);

		// Один корень. f'(x) = 0 или f'(x) > 0
		if ((D_df < -epsilon) || (abs(D_df) <= epsilon))
		{
			roots = new long double[2];
			roots[0] = 1;
			if (abs(c) < epsilon)
			{
				roots[1] = 0;
				return roots;
			}

			else
			{
				if (c > 0)
				{
					roots[1] = bisection(f, INFINITY, 0, epsilon);
					return roots;
				}
				else
				{
					roots[1] = bisection(f, 0, INFINITY, epsilon);
					return roots;
				}
			}

		}

		else
		{
			long double *df_roots = new long double [3];
			df_roots = findRoots(derivative, epsilon);
			long double a = df_roots[1];
			long double b = df_roots[2];


			long double f_a = f->calcF(a);
			long double f_b = f->calcF(b);

			// Один корень на [b; inf)
			if ((f_a < -epsilon) && (f_b < -epsilon))
			{
				roots = new long double[2];
				roots[0] = 1;
				roots[1] = bisection(f, b, INFINITY, epsilon);
				return roots;
			}

			// Один корень на (-inf; a]
			else if ((f_a > epsilon) && (f_b > epsilon))
			{
				roots = new long double[2];
				roots[0] = 1;
				roots[1] = bisection(f, INFINITY, a, epsilon);
				return roots;
			}

			// Три корня
			else if ((f_a > epsilon) && (f_b < -epsilon))
			{
				roots = new long double [4];
				roots[0] = 3;
				roots[1] = bisection(f, INFINITY, a, epsilon);
				roots[2] = bisection(f, a, b, epsilon);
				roots[3] = bisection(f, b, INFINITY, epsilon);
				return roots;
			}

			// Два корня
			else if ((abs(f_b) < epsilon) && (f_a > epsilon))
			{
				roots = new long double [3];
				roots[0] = 2;
				roots[1] = b;
				roots[2] = bisection(f, INFINITY, a, epsilon);
				return roots;
			}

			// Два корня
			else if ((abs(f_a) < epsilon) && (f_b < epsilon))
			{
				roots = new long double [3];
				roots[0] = 2;
				roots[1] = a;
				roots[2] = bisection(f, b, INFINITY, epsilon);
				return roots;
			}

			// Один корень
			else
			{
				roots = new long double [2];
				roots[0] = 1;
				roots[1] = bisection(f, a, b, epsilon);
			}
		}
	}
	return roots;
}

void normaliseInput(long double* coefficients_list, int degree)
{
	if (coefficients_list[0] == 0)
	{
		std::cout << "Коэфициент высшей степени равен нулю" << std::endl;
		exit(-1);
	}
	if (coefficients_list[0] < 0)
	{
		long double d = coefficients_list[0];
		for (int i = 0; i <= degree; i++) {
			coefficients_list[i] = coefficients_list[i] / d;
		}
	}
}

void printRoots(long double *arr)
{
	long double n = arr[0];
	std::cout << "Количество корней = " << n << std::endl << "Корни:" << std::endl;
	for (int i = 1; i <= n; i++)
	{
		std::cout << arr[i] << " ";
	}
	std::cout << std::endl;
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
	normaliseInput(coeficients, n);
	Polynomial* f = new Polynomial(n, coeficients);
	f->printCoefs();
	std::cout << std::endl;
	//long double root_x = bisection(f, -10, 10, 0.05);
	long double *root_x = findRoots(f, 0.1);
	printRoots(root_x);
	delete[] coeficients;
	return 0;
}
