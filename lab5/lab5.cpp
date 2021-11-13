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
		//return pow(3, x) - 2 * x + 5;
}

double diff(double x)
{
	if (variant == 4)
		return 2 * exp(x) - 5;
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
}

int main()
{
	auto difTable = SplitDifTable();

	NewtonMethod(difTable);
	cubicSplain(difTable);

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
