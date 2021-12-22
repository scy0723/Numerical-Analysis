#include <iostream>

#include <math.h>
#define CONV(i) ((float)(i))

using namespace std;

//float machar
void machar(int* ibeta, int* it, int* irnd, int* ngrd, int* machep, int* negep,
	int* iexp, int* minexp, int* maxexp, float* eps, float* epsneg,
	float* xmin, float* xmax)
{
	int i, itemp, iz, j, k, mx, nxres;
	float a, b, beta, betah, betain, one, t, temp, temp1, tempa, two, y, z, zero;

	one = CONV(1);
	two = one + one;
	zero = one - one;
	a = one;
	do {
		a += a;
		temp = a + one;
		temp1 = temp - a;
	} while (temp1 - one == zero);
	b = one;
	do {
		b += b;
		temp = a + b;
		itemp = (int)(temp - a);
	} while (itemp == 0);
	*ibeta = itemp;
	beta = CONV(*ibeta);
	*it = 0;
	b = one;
	do {
		++(*it);
		b *= beta;
		temp = b + one;
		temp1 = temp - b;
	} while (temp1 - one == zero);
	*irnd = 0;
	betah = beta / two;
	temp = a + betah;
	if (temp - a != zero) *irnd = 1;
	tempa = a + beta;
	temp = tempa + betah;
	if (*irnd == 0 && temp - tempa != zero) *irnd = 2;
	*negep = (*it) + 3;
	betain = one / beta;
	a = one;
	for (i = 1; i <= (*negep); i++) a *= betain;
	b = a;
	for (;;) {
		temp = one - a;
		if (temp - one != zero) break;
		a *= beta;
		--(*negep);
	}
	*negep = -(*negep);
	*epsneg = a;
	*machep = -(*it) - 3;
	a = b;
	for (;;) {
		temp = one + a;
		if (temp - one != zero) break;
		a *= beta;
		++(*machep);
	}
	*eps = a;
	*ngrd = 0;
	temp = one + (*eps);
	if (*irnd == 0 && temp * one - one != zero) *ngrd = 1;
	i = 0;
	k = 1;
	z = betain;
	t = one + (*eps);
	nxres = 0;
	for (;;) {
		y = z;
		z = y * y;
		a = z * one;
		temp = z * t;
		if (a + a == zero || fabs(z) >= y) break;
		temp1 = temp * betain;
		if (temp1 * beta == z) break;
		++i;
		k += k;
	}
	if (*ibeta != 10) {
		*iexp = i + 1;
		mx = k + k;
	}
	else {
		*iexp = 2;
		iz = (*ibeta);
		while (k >= iz) {
			iz *= *ibeta;
			++(*iexp);
		}
		mx = iz + iz - 1;
	}
	for (;;) {
		*xmin = y;
		y *= betain;
		a = y * one;
		temp = y * t;
		if (a + a != zero && fabs(y) < *xmin) {
			++k;
			temp1 = temp * betain;
			if (temp1 * beta == y && temp != y) {
				nxres = 3;
				*xmin = y;
				break;
			}
		}
		else break;
	}
	*minexp = -k;
	if (mx <= k + k - 3 && *ibeta != 10) {
		mx += mx;
		++(*iexp);
	}
	*maxexp = mx + (*minexp);
	*irnd += nxres;
	if (*irnd >= 2) *maxexp -= 2;
	i = (*maxexp) + (*minexp);
	if (*ibeta == 2 && !i) --(*maxexp);
	if (i > 20) --(*maxexp);
	if (a != y) *maxexp -= 2;
	*xmax = one - (*epsneg);
	if ((*xmax) * one != *xmax) *xmax = one - beta * (*epsneg);
	*xmax /= (*xmin * beta * beta * beta);
	i = (*maxexp) + (*minexp) + 3;
	for (j = 1; j <= i; j++) {
		if (*ibeta == 2) *xmax += *xmax;
		else *xmax *= beta;
	}
}

//double machar
void machar(int* ibeta, int* it, int* irnd, int* ngrd, int* machep, int* negep,
	int* iexp, int* minexp, int* maxexp, double* eps, double* epsneg,
	double* xmin, double* xmax)
{
	int i, itemp, iz, j, k, mx, nxres;
	double a, b, beta, betah, betain, one, t, temp, temp1, tempa, two, y, z, zero;

	one = CONV(1);
	two = one + one;
	zero = one - one;
	a = one;
	do {
		a += a;
		temp = a + one;
		temp1 = temp - a;
	} while (temp1 - one == zero);
	b = one;
	do {
		b += b;
		temp = a + b;
		itemp = (int)(temp - a);
	} while (itemp == 0);
	*ibeta = itemp;
	beta = CONV(*ibeta);
	*it = 0;
	b = one;
	do {
		++(*it);
		b *= beta;
		temp = b + one;
		temp1 = temp - b;
	} while (temp1 - one == zero);
	*irnd = 0;
	betah = beta / two;
	temp = a + betah;
	if (temp - a != zero) *irnd = 1;
	tempa = a + beta;
	temp = tempa + betah;
	if (*irnd == 0 && temp - tempa != zero) *irnd = 2;
	*negep = (*it) + 3;
	betain = one / beta;
	a = one;
	for (i = 1; i <= (*negep); i++) a *= betain;
	b = a;
	for (;;) {
		temp = one - a;
		if (temp - one != zero) break;
		a *= beta;
		--(*negep);
	}
	*negep = -(*negep);
	*epsneg = a;
	*machep = -(*it) - 3;
	a = b;
	for (;;) {
		temp = one + a;
		if (temp - one != zero) break;
		a *= beta;
		++(*machep);
	}
	*eps = a;
	*ngrd = 0;
	temp = one + (*eps);
	if (*irnd == 0 && temp * one - one != zero) *ngrd = 1;
	i = 0;
	k = 1;
	z = betain;
	t = one + (*eps);
	nxres = 0;
	for (;;) {
		y = z;
		z = y * y;
		a = z * one;
		temp = z * t;
		if (a + a == zero || fabs(z) >= y) break;
		temp1 = temp * betain;
		if (temp1 * beta == z) break;
		++i;
		k += k;
	}
	if (*ibeta != 10) {
		*iexp = i + 1;
		mx = k + k;
	}
	else {
		*iexp = 2;
		iz = (*ibeta);
		while (k >= iz) {
			iz *= *ibeta;
			++(*iexp);
		}
		mx = iz + iz - 1;
	}
	for (;;) {
		*xmin = y;
		y *= betain;
		a = y * one;
		temp = y * t;
		if (a + a != zero && fabs(y) < *xmin) {
			++k;
			temp1 = temp * betain;
			if (temp1 * beta == y && temp != y) {
				nxres = 3;
				*xmin = y;
				break;
			}
		}
		else break;
	}
	*minexp = -k;
	if (mx <= k + k - 3 && *ibeta != 10) {
		mx += mx;
		++(*iexp);
	}
	*maxexp = mx + (*minexp);
	*irnd += nxres;
	if (*irnd >= 2) *maxexp -= 2;
	i = (*maxexp) + (*minexp);
	if (*ibeta == 2 && !i) --(*maxexp);
	if (i > 20) --(*maxexp);
	if (a != y) *maxexp -= 2;
	*xmax = one - (*epsneg);
	if ((*xmax) * one != *xmax) *xmax = one - beta * (*epsneg);
	*xmax /= (*xmin * beta * beta * beta);
	i = (*maxexp) + (*minexp) + 3;
	for (j = 1; j <= i; j++) {
		if (*ibeta == 2) *xmax += *xmax;
		else *xmax *= beta;
	}
}
#undef CONV

//float get_eps 
float get_eps_f() { //1+2^-n=1인 최소의 n을 찾아야 함
	float n = 1.0f; 
	while (1 + n / 2 > 1.0f) //n을 계속 반으로 줄여감
		n /= 2;
	return n ;
}

//double get_eps
double get_eps_d() {//1+2^-n=1인 최소의 n을 찾아야 함
	double n = 1.0f;
	while (1 + n / 2 > 1.0f)//n을 계속 반으로 줄여감
		n /= 2;
	return n ;
}

int main() {
	int ibeta, iexp, irnd, it, machep, maxexp, minexp, negep, ngrd;
	float eps, epsneg, xmax, xmin; //float
	double epsd, epsnegd, xmaxd, xmind; //double 

	machar(&ibeta, &it, &irnd, &ngrd, &machep, &negep, &iexp, &minexp, &maxexp, &eps, &epsneg, &xmin, &xmax);
	cout << "Machar float : " << eps << endl; //float epsilon calculated by machar

	machar(&ibeta, &it, &irnd, &ngrd, &machep, &negep, &iexp, &minexp, &maxexp, &epsd, &epsnegd, &xmind, &xmaxd);
	cout << "Machar double : " << epsd << endl; //double epsilon calculated by machar

	cout << endl;

	cout << "My float eps: " << get_eps_f() << endl; //float epsilon calculated by get_eps_f
	cout << "My double eps : " << get_eps_d() << endl; //double epsilon calculated by get_eps_d

	return 0;
}