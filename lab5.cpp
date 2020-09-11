﻿#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;





void coefC(double* cArr, double** mass, int n, double h)
{
	int i, j;
	double** A = new double* [n - 2];
	for (i = 0; i < n - 2; i++) 
	{
		A[i] = new double[n - 2];
	}
	for (i = 1; i < n - 3; i++) 
	{
		j = i - 1;
		A[i][j] = h;
		A[i][j + 1] = 4 * h;
		A[i][j + 2] = h;
	}
	A[0][0] = 4 * h;
	A[0][1] = h;
	A[n - 3][n - 4] = h;
	A[n - 3][n - 3] = 4 * h;

	double** F = new double* [n - 2];
	for (i = 0; i < n - 2; i++) 
	{
		F[i] = new double[1];
	}

	for (i = 1; i < n - 1; i++) 
	{
		F[i - 1][0] = 3 * ((mass[i + 1][1] - mass[i][1]) / h - (mass[i][1] - mass[i - 1][1]) / h);
	}
	for (i = 1; i < n - 2; i++) 
	{
		F[i][0] -= (A[i][i - 1] * F[i - 1][0] / A[i - 1][i - 1]);
		A[i][i] -= (A[i][i - 1] * A[i - 1][i] / A[i - 1][i - 1]);
	}

	for (i = 1; i < n - 2; i++) 
	{
		F[i][0] -= (A[i][i - 1] * F[i - 1][0] / A[i - 1][i - 1]);
		A[i][i] -= (A[i][i - 1] * A[i - 1][i] / A[i - 1][i - 1]);
	}

	cArr[n - 2] = F[n - 3][0] / A[n - 3][n - 3];
	for (i = n - 4; i >= 0; i--) 
	{
		cArr[i + 1] = (F[i][0] - A[i][i + 1] * cArr[i + 2]) / A[i][i];
	}
	cArr[0] = 0;
	cArr[n - 1] = 0;
	for (int i = 0; i < n - 2; i++) 
	{
		delete[] A[i];
		delete[] F[i];
	}
	delete[] A;
	delete[] F;
}



double splineFunc(double x, double* a, double* b, double* c, double* d, double** func, int n)
{
	int i = 0;
	for (int j = 1; j < n - 1; j++)
	{
		if (x > func[j][0])
			i++;
		else
			break;
	}
	return a[i] + b[i] * (x - func[i][0]) + c[i] * pow((x - func[i][0]), 2) + d[i] * pow((x - func[i][0]), 3);
}


int main(void)
{
	using namespace std;
	int i, n = 9;
	double** func;
	func = new double* [n];
	double* a = new double[n - 1];
	double* b = new double[n - 1];
	double* c = new double[n];
	double* d = new double[n - 1];

	for (int i = 0; i <= n - 1; i++)
	{
		func[i] = new double[n];
		for (int j = 0; j <= 1; j++)
		{
			if ((i == 0) && (j == 0))
				func[i][j] = -1;
			else if ((i == 0) && (j == 1))
				func[i][j] = -0.333;
			else if ((i == 1) && (j == 0))
				func[i][j] = -0.75;
			else if ((i == 1) && (j == 1))
				func[i][j] = 0;
			else if ((i == 2) && (j == 0))
				func[i][j] = -0.5;
			else if ((i == 2) && (j == 1))
				func[i][j] = -0.125;
			else if ((i == 3) && (j == 0))
				func[i][j] = -0.25;
			else if ((i == 3) && (j == 1))
				func[i][j] = -0.056;
			else if ((i == 4) && (j == 0))
				func[i][j] = 0;
			else if ((i == 4) && (j == 1))
				func[i][j] = 0;
			else if ((i == 5) && (j == 0))
				func[i][j] = 0.25;
			else if ((i == 5) && (j == 1))
				func[i][j] = 0.046;
			else if ((i == 6) && (j == 0))
				func[i][j] = 0.5;
			else if ((i == 6) && (j == 1))
				func[i][j] = 0.083;
			else if ((i == 7) && (j == 0))
				func[i][j] = 0.75;
			else if ((i == 7) && (j == 1))
				func[i][j] = 0.115;
			else if ((i == 8) && (j == 0))
				func[i][j] = 1;
			else if ((i == 8) && (j == 1))
				func[i][j] = 0.143;
		
		}
	}
	double h = 0.125;
	coefC(c, func, n, h);
	for (i = 0; i < n - 1; i++)
	{
		a[i] = func[i][1];
		b[i] = (func[i + 1][1] - func[i][1]) / h - c[i] * h - (c[i + 1] - c[i]) * h / 3.;
		d[i] = (c[i + 1] - c[i]) / (3 * h);
	}




	double up = 0, down = 1;
	int ndrugaai = 999999;
	double width =(down - up) / ndrugaai;


	double trapezoidal_integral = 0;
	for (int step = 0; step < ndrugaai; step++) {
		const double x1 = up + step * width;
		const double x2 = up + (step + 1) * width;

		trapezoidal_integral += 0.5 * (x2 - x1) * (splineFunc(x1, a, b, c, d, func, n) + splineFunc(x2, a, b, c, d, func, n));
	}

	double trapeciya;
	double tochnost_trapecii;
	double polovina_trapeciya;
	double ndrugaaitrap;
		trapeciya = trapezoidal_integral;
		ndrugaai = 999999;
		ndrugaaitrap = ndrugaai / 2;
		width = (down - up) / ndrugaaitrap;
		trapezoidal_integral = 0;
		for (int step = 0; step < ndrugaaitrap; step++)
		{
			const double x1 = up + step * width;
			const double x2 = up + (step + 1) * width;
			trapezoidal_integral += 0.5 * (x2 - x1) * (splineFunc(x1, a, b, c, d, func, n) + splineFunc(x2, a, b, c, d, func, n));
		}
		polovina_trapeciya = trapezoidal_integral;
		tochnost_trapecii = fabs((trapeciya - polovina_trapeciya) / trapeciya);

	cout << "Trapezoid method: " << trapezoidal_integral << "\t" << "Runge's rule for trapezoid method: " << tochnost_trapecii << "\t" << endl;

	ofstream trap("ans1.dat");
	trap << trapeciya << "\t" << tochnost_trapecii << endl;
	trap.close();






	return 0;
}
