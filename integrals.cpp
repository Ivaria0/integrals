#include <iostream>
#include <math.h>

using namespace std;

double a = -4;
double b = 4;

double f(double x)
{
	double r;
	r = 1 / (1 + exp(x));
	return r;
}

double Simpson(int n)
{
	double res = 0, r1, r2 = a, f1, f2 = f(r2);
	for (int i = 1; i <= n; i++) {
		r1 = r2;
		r2 += (b - a) / n;
		f1 = f2;
		f2 = f(r2);
		res += f1 + 4 * f((r2 + r1) / 2) + f2;
	}
	res *= ((b - a) / (6 * n));
	return res;
}

double Gauss2(int n)
{
	double res = 0, r1, r2 = a, x1, x2;
	for (int i = 1; i <= n; i++) {
		r1 = r2;
		r2 += (b - a) / n;
		x1 = (r1 + r2) / 2 - sqrt(3) * (r2 - r1) / 6;
		x2 = (r1 + r2) / 2 + sqrt(3) * (r2 - r1) / 6;
		res += f(x1) + f(x2);
	}
	res *= ((b - a) / (2 * n));
	return res;
}

double integral()
{
	double res;
	res = log(abs(exp(b) / (exp(b) + 1))) - log(abs(exp(a) / (exp(a) + 1)));
	return res;
}

double RungeS(int n, int p)
{
	double res;
	res = abs(Simpson(n) - Simpson(n / 2)) / (pow(2, p) - 1);
	return res;
}

double RungeG(int n, int p)
{
	double res;
	res = abs(Gauss2(n) - Gauss2(n / 2)) / (pow(2, p) - 1);
	return res;
}

void main()
{
	int n;
	cin >> n;
	cout << "Simpson(" << n << ")=" << Simpson(n) << endl;
	cout << "Gauss(" << n << ")=" << Gauss2(n) << endl;
  
	double maxS = 0, maxG = 0;
	int nS = 0, nG = 0;
  
	for (int i = 1; i <= 20; i++) {
		n = pow(2, i);
		cout << "Simpson(" << n << ")=" << Simpson(n) << endl;
		cout << "Gauss(" << n << ")=" << Gauss2(n) << endl;
		cout << "Runge for Simpson(" << n << ")=" << RungeS(n, 3) << endl;
		cout << "Runge for Gauss(" << n << ")=" << RungeG(n, 3) << endl;
    
		if (RungeS(n, 3) > maxS) {
			maxS = RungeS(n, 3);
			nS = n;
		}
    
		if (RungeG(n, 3) > maxG) {
			maxG = RungeG(n, 3);
			nG = n;
		}
	}
  
	int pS = 1, pG = 1;
  
	while (abs(integral() - Simpson(n)) < RungeS(n, pS)) pS++;
  
	while (abs(integral() - Gauss2(n)) < RungeG(n, pG)) pG++;
  
	cout << "Order of accuracy of Simpson: " << pS << endl;
	cout << "Order of accuracy of Gauss: " << pG << endl;
}
