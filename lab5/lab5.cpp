// lab5.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <vector>
#include <iomanip>
#include "MatrixOperation.h"
#include "VectorsOperation.h"

using namespace std;

int variant = 4;
int x0 = 1;
int xn = 2;
int n = 5;
double h = (double)(xn - x0) / n;

double f(double x)
{
	if (variant == 4)
		return 2 * exp(x) - 5 * x;
	return exp(x) + x + 1;
}

double diff(double x)
{
	if (variant == 4)
		return 2 * exp(x) - 5;
	return exp(x) + 1;
}


double phi0(double tau)
{
	return (1 + 2 * tau) * pow(1 - tau, 2);
}

double phi1(double tau)
{
	return tau * pow(1 - tau, 2);
}

vector<double> sqrtMethod(vector<vector<double>> a, vector<double> b)
{
	vector<vector<double>> S(a.size(), (vector<double>(a.size())));
	vector<vector<double>> d(a.size(), (vector<double>(b.size())));


	for (int i = 0; i < S.size(); ++i)
	{
		double sum = 0;

		for (int k = 0; k < i; ++k)
			sum += S[k][i] * S[k][i] * d[k][k];

		d[i][i] = copysign(1.0, a[i][i] - sum);
		S[i][i] = sqrt(a[i][i] - sum);

		double temp = 1 / (S[i][i] * d[i][i]);

		for (int j = i + 1; j < S.size(); ++j)
		{
			sum = 0;

			for (int k = 0; k < i; ++k)
				sum += S[k][i] * S[k][j] * d[k][k];

			S[i][j] = (a[i][j] - sum) * temp;
		}

	}

	vector<double> y(b.size());
	vector<double> x(b.size());

	y[0] = b[0] / (S[0][0] * d[0][0]);

	for (int i = 1; i < b.size(); ++i)
	{
		double sum = 0;
		for (int j = 0; j < i; ++j)
			sum += S[j][i] * y[j] * d[j][j];

		y[i] = (b[i] - sum) / (S[i][i] * d[i][i]);
	}

	x[b.size() - 1] = y[b.size() - 1] / S[b.size() - 1][b.size() - 1];

	for (int i = x.size() - 2; i > -1; --i)
	{
		double sum = 0;
		for (int j = i; j < x.size(); ++j)
			sum += S[i][j] * x[j];

		x[i] = (y[i] - sum) / S[i][i];
	}

	return x;
}

vector<vector<double>> SplitDifTable()
{
	vector<vector<double>> res(n + 1, (vector<double>(n + 2)));

	for(int i=0;i<res.size();++i)
	{
		double x = x0 + i * h;
		res[i][0] = x;
		res[i][1] = f(x);
	}

	for(int i=2;i<=res.size();++i)
		for(int j=0;j<=res.size()-i;++j)
		{
			res[j][i] = (res[j + 1][i - 1] - res[j][i - 1]) / (res[i + j - 1][0] - res[j][0]);
		}

	return res;
}


vector<vector<double>> SplitDifTable2()
{
	vector<vector<double>> res(n + 1, (vector<double>(n + 2)));

	for (int i = 0; i < res.size(); ++i)
	{
		double x = x0 + i * h;
		res[i][0] = f(x);
		res[i][1] = x;
	}

	for (int i = 2; i <= res.size(); ++i)
		for (int j = 0; j <= res.size() - i; ++j)
		{
			res[j][i] = (res[j + 1][i - 1] - res[j][i - 1]) / (res[i + j - 1][0] - res[j][0]);
		}

	return res;
}

void NewtonMethod(vector<vector<double>> difTable)
{
	int fact = 1;

	for (int i = 2; i < n + 2; ++i)
		fact *= i;

	double M6;

	if (variant == 4)
		M6 = 2 * exp(2);
		//M6 = 9 * pow(log(3), 6);

	cout << "M6=" << fixed << setprecision(7) <<M6<< '\n';

	for(int i=0;i<n;++i)
	{
		double x = x0 + (i + 0.5) * h;
		double Pn = difTable[0][1];
		double w = 1;

		for(int j=1;j<n+1;++j)
		{
			w *= x - difTable[j-1][0];
			Pn += w * difTable[0][j+1];
		}

		w *= x - difTable[5][0];

		cout << fixed << setprecision(2) << x << " " << setprecision(5) << f(x) << " " << Pn << " "<<scientific<<abs(Pn-f(x))<<" "<<M6*abs(w)/fact<<'\n';

	}
}

void cubicSplain(vector<vector<double>> difTable)
{
	vector<double> m(n+1);

	m[0] = diff(1);
	m[n] = diff(2);

	vector<double> alpha(n+1);
	vector<double> beta(n+1);

	alpha[1] = 0;
	beta[1] = diff(1);

	for(int i=1;i<n;++i)
	{
		alpha[i + 1] = -1 / (4 + alpha[i]);
		beta[i + 1] = (3 * (difTable[i + 1][1] - difTable[i - 1][1]) / h - beta[i]) / (4 + alpha[i]);
	}

	for (int i = n - 1; i > -1; --i)
		m[i] = alpha[i + 1] * m[i + 1] + beta[i + 1];

	double M4, M5;

	if(variant==4)
	{
		M4 = M5 = 2 * exp(2);
	}

	cout << "M5=" << fixed << setprecision(7) << M5 << '\n';
	cout << "M4=" << fixed << setprecision(7) << M4 << '\n';

	for(int i=0;i<n+1;++i)
	{
		double x = x0 + i * h;

		cout << fixed << setprecision(1) << x << " " << setprecision(6) << diff(x) << " " << m[i] << " " << abs(diff(x) - m[i]) << " " << M5 / 60 * pow(h, 4) << endl;
	}

	double tau = 0.5;//(x-xi)/h const in the fact

	for(int i=0;i<n;++i)
	{
		double x = x0 + (i + 0.5) * h;
		double s = phi0(tau) * difTable[i][1] + phi0(1 - tau) * difTable[i+1][1] + h * (phi1(tau) * m[i] - phi1(1 - tau) * m[i + 1]);

		cout << fixed << setprecision(2) << x << " " << setprecision(6) << f(x) << " " << s << " " << abs(f(x) - s) << " " << (M4/384+M5*h/240)*pow(h,4) << '\n';
	}
}

void DiscreteRMS(vector<vector<double>> difTable)
{
	vector<vector<double>> a(3, vector<double>(3));

	for(int i=0;i<n;++i)
	{
		a[0][0]++;
		a[0][1] += difTable[i][0];
		a[0][2] += pow(difTable[i][0], 2);
		a[1][2] += pow(difTable[i][0], 3);
		a[2][2] += pow(difTable[i][0], 4);
	}

	a[1][0] = a[0][1];
	a[1][1] = a[0][2];
	a[2][0] = a[0][2];
	a[2][1] = a[1][2];

	cout << "Матрица\n";

	for(int i=0;i<a.size();++i)
	{
		for (int j = 0; j < a[i].size(); ++j)
			cout << fixed << setprecision(6) << a[i][j] << " ";
		cout << '\n';
	}

	vector<double> b(3);

	for(int i=0;i<n;++i)
	{
		b[0] += f(difTable[i][0]);
		b[1] += f(difTable[i][0]) * difTable[i][0];
		b[2] += f(difTable[i][0]) * difTable[i][0] * difTable[i][0];
	}

	cout << "Вектор правых частей\n";

	for (int i = 0; i < b.size(); ++i)
		cout << fixed << setprecision(6) << b[i] << " ";

	cout << '\n';

	vector<double> solvedX = sqrtMethod(a, b);

	//У меня сошлось с графиками
	cout << fixed << setprecision(5) << "P2(x)= " << solvedX[0] << " +(" << solvedX[1] << "*x) + (" << solvedX[2] << "*x^2)\n";

	double normF = 0;
	double normG = 0;

	for(int i=0;i<n;++i)
	{
		normF += pow(difTable[i][1], 2);
		normG += pow(solvedX[0] + solvedX[1] * difTable[i][0] + solvedX[2] * pow(difTable[i][0], 2), 2);
	}

	cout << "Error= " << fixed << scientific << sqrt(abs(normF - normG))<<'\n';

	for(int i=0;i<n;++i)
	{
		double x = x0 + i * h;

		cout <<fixed<<setprecision(2)<< x << " " << fixed << setprecision(15) << f(x) - (solvedX[0] + solvedX[1] * x + solvedX[2] * x * x) << '\n';
	}
}

void IntegretedRMS(vector<vector<double>> difTable)
{
	vector<vector<double>> a(3, (vector<double>(3)));

	/*Тут хитрый алгоритм какой - то, чуть раздуплил
	 * У нас есть полином второй степени, то есть 1 x x^2
	 * Считаем интегралы от 1 до 2 всех парных комбинаций типа g1*g1 g1*g2 g1*g3 g2*g1 g2*g2 ...
	 * И соответсвенно результат сюда вносим
	 */
	a[0][0] = 1;
	a[0][1] = a[1][0] = 3 / 2.;
	a[0][2] = a[1][1] = a[0][2] = 7 / 3.;
	a[1][2]=a[2][1] = 15 / 4.;
	a[2][2] = 31 / 5.;

	for (int i = 0; i < a.size(); ++i)
	{
		for (int j = 0; j < a[i].size(); ++j)
			cout << fixed << setprecision(6) << a[i][j] << " ";
		cout << '\n';
	}

	vector<double> b(3);

	/* Тут тоже самое почти, только мы умножаем нашу f на g => f*g1, f*g2, f*g3
	 * От этого всего интеграл берем от 1 до 2
	 * Значения готовы
	 */

	if(variant==4)
	{
		b[0] = 2 * exp(2)- 2 * exp(1) - 15 / 2.;
		b[1] = 2 * exp(2) - 35 / 3.;
		b[2] = 4 * exp(2) - 2 * exp(1) - 75 / 4.;
	}

	for (int i = 0; i < b.size(); ++i)
		cout << fixed << setprecision(6) << b[i] << " ";

	cout << '\n';

	vector<double> solvedX = sqrtMethod(a, b);

	//У меня сошлось с графиками
	cout << fixed << setprecision(5) << "P2(x)= " << solvedX[0] << " +(" << solvedX[1] << "*x) + (" << solvedX[2] << "*x^2)\n";

	/*
	 * Тут руками надо считать int(f^2)-int(P2(x)^2)
	 * а после корень из этого взять
	 */
	double normF;
	double normG;

	if (variant == 4)
	{
		normF = (6*exp(4)-66*exp(2)+175)/3.;
		normG = 745438076087 /150000000000.;
	}

	cout << "Error= " << fixed << scientific << sqrt(abs(normF - normG)) << '\n';

	for (int i = 0; i < n; ++i)
	{
		double x = x0 + i * h;

		cout << fixed << setprecision(2) << x << " " << fixed << setprecision(15) << f(x) - (solvedX[0] + solvedX[1] * x + solvedX[2] * x * x) << '\n';
	}
}


double FindRoot(double A) 
{
	double a = x0;
	double b = xn;
	double c;
	while (abs((b - a) / 2) > eps) {
		c = (a + b) / 2;
		if (((diff(a) - A) * (diff(c) - A)) > 0) a = c;
		else b = c;
	}
	return c;
}


void Poli_one()
{
	double a0, a1, d;
	a1 = (f(xn) - f(x0)) / (xn - x0);
	d = FindRoot(a1);
	a0 = (f(x0) + f(d) - a1 * (x0 + d)) / 2;

	cout << "P1(x) = " << fixed << setprecision(4) << a0 << " + " << a1 << " * x,		d = "<<d<<endl;
	cout << "L(a) = " << f(x0) - (a0 + a1 * x0) << endl;
	cout << "L(d) = " << f(d) - (a0 + a1 * d) << endl;
	cout << "L(b) = " << f(xn) - (a0 + a1 * xn) << endl;

	cout << setw(4) << 'X' <<setw(19)<< "Погрешность" << endl;
	for (int i = 0; i < 6; i++)
	{
		cout <<setw(4)<<setprecision(1)<< x0 + i * 0.2 <<"     "<<setprecision(15)<<setw(17)<< f(x0 + i * 0.2) - (a0 + a1 * (x0 + i * 0.2)) << endl;
	}
}

void opposite_RMS()
{
	vector<vector<double>> dif_matr = SplitDifTable2();

	for (int i = 0; i < dif_matr.size(); i++)
	{
		for (int j = 0; j <dif_matr[i].size(); j++)
		{
			cout <<fixed<<setprecision(5)<<setw(10)<< dif_matr[i][j] << " ";
		}
		cout << endl;
	}
	double C = (f(x0)+f(xn))/2;
	cout << "C = " << C<<endl;


	double w = 1;
	double answ = dif_matr[0][1];
	for (int i = 1; i < n + 1; i++)
	{
		w *= (C - dif_matr[i-1][0]);
		answ += w * dif_matr[0][i+1];
	}

	cout <<setprecision(15)<< "Корень: " << answ << endl;
	cout << "Невязка " << scientific<< abs(f(answ) - C);
}

int main()
{
	auto difTable = SplitDifTable();

	NewtonMethod(difTable);
	cubicSplain(difTable);
	DiscreteRMS(difTable);
	IntegretedRMS(difTable);
	Poli_one();
	opposite_RMS();

}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
