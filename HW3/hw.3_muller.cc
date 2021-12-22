#include<stdio.h>
#include<math.h>
#include<iostream>
using namespace std;

const int MAXIT= 10000;

float f(float x)
{
	return pow(x, 3) + 4 * x * x + 15 * x - 27;
}

void Muller(float a, float b, float c)
{
	int i;
	float res;

	for (i = 0;; ++i)
	{

		float f1 = f(a);
		float f2 = f(b);
		float f3 = f(c);
		float d1 = f1 - f3;
		float d2 = f2 - f3;
		float h1 = a - c;
		float h2 = b - c;
		float a0 = f3;
		float a1 = (((d2 * pow(h1, 2)) - (d1 * pow(h2, 2))) / ((h1 * h2) * (h1 - h2)));
		float a2 = (((d1 * h2) - (d2 * h1)) / ((h1 * h2) * (h1 - h2)));
		float x = ((-2 * a0) / (a1 + abs(sqrt(a1 * a1 - 4 * a0 * a2))));
		float y = ((-2 * a0) / (a1 - abs(sqrt(a1 * a1 - 4 * a0 * a2))));

		if (x >= y)
			res = x + c;
		else
			res = y + c;

		float m = res * 100;
		float n = c * 100;
		m = floor(m);
		n = floor(n);
		if (m == n)
			break;
		a = b;
		b = c;
		c = res;
		if (i > MAXIT)
		{
			break;
		}
	}
	cout << "root of function:  " << res;
}

// Driver main function
int main()
{
	printf("Muller method: ex) x^3 + 4x^2 + 15x - 27\n");
	float a = 0, b = 1, c = 2;
	Muller(a, b, c);
	return 0;
}
