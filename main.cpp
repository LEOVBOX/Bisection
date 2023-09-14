#include <iostream>
#include <stdarg.h>
#include <math.h>
#include <vector>

class Root
{
public:
	long double x;
	int multiplicity;
	Root(long double x, int multiplicity)
	{
		this->x = x;
		this->multiplicity = multiplicity;
	}

	Root(long double x)
	{
		this->x = x;
		this->multiplicity = 1;
	}

	Root()
	{
		this->x = 0;
		this->multiplicity = 0;
	}
};

class Polynomial
{
	// p_0 * x^n + p_1 * x^(n-1) + ... p_(n-1)*x + p_n
	// p_0 > 0
public:
	long double* coefficients_list;
	int degree;

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
	while (f_x > 0)
	{
		x_L -= delta;
		f_x = f->calcF(x_L);
	}
	return x_L;
}

long double findXR(Polynomial* f, long double x_L)
{
	long double delta = 1;
	long double x_R = x_L + delta;
	long double f_x = f->calcF(x_R);
	while (f_x < 0)
	{
		x_R += delta;
		f_x = f->calcF(x_R);
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
	Polynomial *df = f->calcDerivative();
	long double df_x = df->calcF(x_average);
	while (abs(f_x) >= epsilon)
	{
		// Ищем решение на [x_average; x_R]
		if ((f_x < -epsilon && df_x > 0) || (f_x > epsilon && df_x < 0))
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
		df_x = df->calcF(x_average);
	}
	return x_average;
}

long double calcD(long double a, long double b, long double c)
{
	return pow(b, 2) - 4 * a * c;
}

std::vector<Root*> findRoots(Polynomial* f, long double epsilon)
{
	std::vector<Root*> solution;

	// Полином первой степени
	if (f->degree == 1)
	{
		long double k = f->coefficients_list[0];
		long double b = f->coefficients_list[1];
		if (k == 0)
		{
			std::cout << "на вход подана константа";
			exit(-1);
		}
		else
		{
			Root *root = new Root(-(b/k));
			solution.push_back(root);
			return solution;
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
			Root *root = new Root(-(b/(2*a)), 2);
			solution.push_back(root);
			return solution;
		}
		else if (D > 0)
		{

			Root* x_1 = new Root((-b - sqrt(D))/(2*a));
			solution.push_back(x_1);
			Root* x_2  = new Root((-b + sqrt(D))/(2*a));
			solution.push_back(x_2);
			return solution;
		}
		else
		{
			std::cout << "уравнение не имеет действительных корней";
			exit(-1);
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
			if (abs(c) < epsilon)
			{
				Root *root = new Root(0);
				solution.push_back(root);
				return solution;
			}

			long double a = f->coefficients_list[1];
			long double b = f->coefficients_list[2];

			// x^3 + 3*z*x^2 + 3*x*z^2 + z^3 = (x + z)^3
			long double z = cbrt(c);
			if (abs(a)/3 == abs(z) && sqrt(abs(b)/3) == abs(z))
			{
				Root* root  = new Root(-z, 3);
				solution.push_back(root);
				return solution;
			}

			else
			{
				if (c > 0)
				{
					Root* root  = new Root(bisection(f, INFINITY, 0, epsilon));
					solution.push_back(root);
					return solution;
				}
				else
				{
					Root* root = new Root(bisection(f, 0, INFINITY, epsilon));
					solution.push_back(root);
					return solution;
				}
			}

		}

		else
		{
			std::vector<Root*> df_solution;
			df_solution = findRoots(derivative, epsilon);
			long double a = df_solution[0]->x;
			long double b = df_solution[1]->x;


			long double f_a = f->calcF(a);
			long double f_b = f->calcF(b);

			// Один корень на [b; inf)
			if ((f_a < -epsilon) && (f_b < -epsilon))
			{
				Root* root = new Root(bisection(f, b, INFINITY, epsilon));
				solution.push_back(root);
				return solution;
			}

			// Один корень на (-inf; a]
			else if ((f_a > epsilon) && (f_b > epsilon))
			{
				Root *root = new Root(bisection(f, INFINITY, a, epsilon));
				solution.push_back(root);
				return solution;
			}

			// Три корня
			else if ((f_a > epsilon) && (f_b < -epsilon))
			{
				Root* root_1 = new Root(bisection(f, INFINITY, a, epsilon));
				Root* root_2 = new Root(bisection(f, a, b, epsilon));
				Root* root_3 = new Root(bisection(f, b, INFINITY, epsilon));
				solution.push_back(root_1);
				solution.push_back(root_2);
				solution.push_back(root_3);
				return solution;
			}

			// x^3 + 3*z*x^2 + 3*x*z^2 + z^3 = (x + z)^3

			// Два корня
			else if ((abs(f_b) < epsilon) && (f_a > epsilon))
			{
				Root *root_1 = new Root(b, 2);
				solution.push_back(root_1);
				Root *root_2 = new Root(bisection(f, INFINITY, a, epsilon));
				solution.push_back(root_2);
				return solution;
			}

			// x^3 + 3*z*x^2 + 3*x*z^2 + z^3 = (x + z)^3

			// Два корня
			else if ((abs(f_a) < epsilon) && (f_b < epsilon))
			{
				Root *root_1 = new Root(a, 2);
				solution.push_back(root_1);
				Root *root_2 = new Root(bisection(f, b, INFINITY, epsilon));
				solution.push_back(root_2);
				return solution;
			}

			// Один корень
			else if (abs(f_a) < epsilon && abs(f_b) < epsilon)
			{
				Root *root = new Root(bisection(f, a, b, epsilon), 3);
				solution.push_back(root);
				return solution;
			}
		}
	}

	return solution;
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

void printRoots(std::vector<Root*> &roots)
{
	size_t n = roots.size();
	std::cout << "Количество корней = " << n << std::endl << "Корни:" << std::endl;
	for (int i = 0; i < n; i++)
	{
		std::cout << "x_" << i << " = " << roots[i]->x << " Кратность = " << roots[i]->multiplicity << std::endl;
	}
	std::cout << std::endl;
}

long double* getCoefficients()
{
	long double* coefficients_list = new long double[4];
	std::cout << "f(x) = x^3 + ax^2 + bx + c" << std::endl;
	coefficients_list[0] = 1;
	std::cout << "a = ";
	std::cin >> coefficients_list[1];
	std::cout << "b = ";
	std::cin >> coefficients_list[2];
	std::cout << "c = ";
	std::cin >> coefficients_list[3];
	std::cout << "Решение уравнениия f(x) = 0" << std::endl;
	return coefficients_list;
}


int main()
{
	int n = 3;
	long double epsilon;
	std::cout << "epsilon = ";
	std::cin >> epsilon;

	auto *coeficients = getCoefficients();
	//normaliseInput(coeficients, n);
	Polynomial* f = new Polynomial(n, coeficients);
	std::cout << std::endl;

	std::vector<Root*> solution = findRoots(f, epsilon);
	printRoots(solution);

	return 0;
}
