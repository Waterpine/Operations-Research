#include "dualsimplex.h"
void DualSimplexMatrix::dual()
{
	int temp;
	double *temp_double;

	temp = this->m;
	this->m = this->n;
	this->n = temp;

	delete[] d;
	delete[] e;
	this->e = new int[this->n];
	this->d = new int[this->m];
	for (int i = 0; i < this->n; i++)
	{
		this->e[i] = 0;
	}
	for (int i = 0; i < this->m; i++)
	{
		this->d[i] = 1;
	}

	temp_double = this->c;
	this->c = this->b;
	this->b = temp_double;

	for (int i = 0; i < this->n; i++)
	{
		this->c[i] = -this->c[i];
	}
	this->dual_a = new double *[this->m];
	for (int i = 0; i < this->m; i++)
	{
		this->dual_a[i] = new double[this->n];
		for (int j = 0; j < this->n; j++)
		{
			this->dual_a[i][j] = this->a[j][i];
		}
	}
	this->a = this->dual_a;
}

void DualSimplexMatrix::prepare()
{
	this->z = 0;
	//新的矩阵
	vector<double> aa[MAX_M];
	//矩阵初始化
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			aa[i].push_back(a[i][j]);
	//处理e
	for (int i = 0; i < n; i++)
	{
		//xi <= 0 变 xi >= 0 改变a和c
		if (e[i] < 0)
		{
			for (int j = 0; j < m; j++)
				aa[j][i] = -aa[j][i];
			c[i] = -c[i];
		}
		//xi无限制变为xi - xk  xi >= 0 && xk >= 0 
		else if (e[i] == 0)
			for (int j = 0; j < m; j++)
				aa[j].push_back(-1);
	}

	//处理d，添加松弛变量
	for (int i = 0; i < m; i++)
	{
		if (d[i] > 0)
		{
			for (int j = 0; j < m; j++)
			{
				if (j == i)
					aa[j].push_back(-1);
				else
					aa[j].push_back(0);
			}
			for (int j = 0; j < aa[i].size(); j++)
			{
				aa[i][j] = -aa[i][j];
			}
			b[i] = -b[i];
		}
		else if (d[i] < 0)
		{
			for (int j = 0; j < m; j++)
			{
				if (j == i)
					aa[j].push_back(1);
				else
					aa[j].push_back(0);
			}
		}
	}

	//a
	int new_n = aa[0].size();
	for (int i = 0; i < m; i++)
	{
		delete[] this->a[i];
		this->a[i] = new double[new_n];
		for (int j = 0; j < new_n; j++)
		{
			a[i][j] = aa[i][j];
		}
	}
	//c
	double *tmp = new double[new_n];
	for (int i = 0; i < new_n; i++)
	{
		if (i < this->n)
			tmp[i] = this->c[i];
		else
			tmp[i] = 0;
	}
	delete[]this->c;
	this->c = tmp;

	//add x
	x = new double[new_n];
	for (int i = 0; i < new_n; i++)
		x[i] = 0;
	this->n = new_n;
}

void DualSimplexMatrix::readFile(const char* path)
{
	ifstream file(path);
	if (!file.is_open())
	{
		cout << "Open file error" << endl;
	}
	file >> this->n >> this->m;
	this->a = new double*[this->m];
	for (int i = 0; i < this->m; i++)
	{
		this->a[i] = new double[this->n];
	}
	this->b = new double[this->m];
	this->c = new double[this->n];
	this->d = new int[this->m];
	this->e = new int[this->n];
	for (int i = 0; i < this->n; i++)
	{
		file >> this->c[i];
	}
	for (int i = 0; i < this->m; i++)
	{
		for (int j = 0; j < this->n; j++)
		{
			file >> this->a[i][j];
		}
		file >> this->b[i];
		file >> this->d[i];
	}
	for (int i = 0; i < this->n; i++)
	{
		file >> this->e[i];
	}
}

void DualSimplexMatrix::print()
{
	cout << "m = " << m << endl;
	cout << "n = " << n << endl;
	cout << "c = ";
	for (int i = 0; i < n; i++)
	{
		cout << c[i] << " ";
	}
	cout << endl;
	cout << "Matrix: " << endl;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << a[i][j] << "\t";
		}
		cout << b[i];
		cout << endl;
	}
}

void DualSimplexMatrix::print_root()
{
	if (x == NULL)
	{
		cout << "No root" << endl;
		return;
	}
	//cout << "The root is: ";
	for (int i = 0; i < this->n; i++)
	{
		cout << x[i] << " ";
	}
	cout << endl;
}

int DualSimplexMatrix::do_counting()
{
	double **abc;
	double **n_abc;
	abc = new double*[m + 1];
	n_abc = new double*[m + 1];
	for (int i = 0; i < m + 1; i++) {
		abc[i] = new double[n];
		n_abc[i] = new double[n];
	}

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			abc[i][j] = a[i][j];
		}
	}
	for (int i = 0; i < m; i++) {
		abc[i][n] = b[i];
	}
	for (int i = 0; i < n; i++) {
		abc[m][i] = c[i];
	}
	abc[m][n] = -z;

	int in_base, out_base;
	int argmin_outbase;
	double min_outbase = 0;
	double temp;
	double cur_min;
	bool is_base;
	while (true) {
		out_base = -1;
		cur_min = 0;
		for (int i = 0; i < m; i++) {
			if (abc[i][n] < cur_min) {
				out_base = i;
				cur_min = abc[i][n];
			}
		}
		if (out_base < 0) {
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < n; j++) {
					if (abc[i][j] == 1) {
						is_base = true;
						for (int k = 0; k < m; k++) {
							if (k != i&&abc[k][j] != 0) {
								is_base = false;
								break;
							}
						}
						if (is_base == true) {
							x[j] = abc[i][n];
						}
					}
				}
			}
			z = -abc[m][n];

			for (int i = 0; i < m; i++) {
				for (int j = 0; j < n; j++) {
					a[i][j] = abc[i][j];
				}
			}

			for (int i = 0; i < m; i++) {
				b[i] = abc[i][n];
			}
			return 1;
		}
		else {
			in_base = -1;
			cur_min = DBL_MAX;
			temp = DBL_MAX;
			for (int i = 0; i < n; i++) {
				if (abc[out_base][i] < 0) {
					temp = fabs(abc[m][i] / abc[out_base][i]);
				}
				if (temp < cur_min) {
					cur_min = temp;
					in_base = i;
				}
			}
			if (in_base < 0) {
				return -1;
			}
			else {
				for (int i = 0; i < m + 1; i++) {
					for (int j = 0; j < n + 1; j++) {
						if (i != out_base) {
							n_abc[i][j] = abc[i][j] - (abc[i][in_base] / abc[out_base][in_base])*abc[out_base][j];
						}
						else {
							n_abc[i][j] = abc[i][j] / abc[out_base][in_base];
						}
					}
				}
				for (int i = 0; i < m + 1; i++) {
					for (int j = 0; j < n + 1; j++) {
						abc[i][j] = n_abc[i][j];
						//cout << abc[i][j] << '\t';
					}
					//cout << endl;
				}
				//cout << endl;
			}
		}
	}
}

void DualSimplexMatrix::make_normal()
{
	for (int i = 0; i < this->m; i++)
	{
		for (int j = 0; j < this->n; j++)
		{
			a[i][j] = -a[i][j];
		}
		b[i] = -b[i];
	}
}
