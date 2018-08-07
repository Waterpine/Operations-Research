#include "simplex.h"

void SimplexMatrix::prepare()
{
	//对<0的b进行变换
	for (int i = 0; i < m; i++)
	{
		if (b[i] < 0)
		{
			b[i] = -b[i];
			d[i] = -d[i];
			for (int j = 0; j < n; j++)
			{
				a[i][j] = -a[i][j];
			}
		}
	}
	//新的矩阵
	vector<double> aa[MAX_M];
	vector<double> aa1[MAX_M];
	//矩阵初始化
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			aa[i].push_back(a[i][j]);
		}
	}
	//处理e
	for (int i = 0; i < n; i++)
	{
		//xi <= 0 变 xi >= 0 改变a和c
		if (e[i] < 0)
		{
			for (int j = 0; j < m; j++)
			{
				aa[j][i] = -aa[j][i];
			}
			c[i] = -c[i];
		}
		//xi无限制变为xi - xk  xi >= 0 && xk >= 0 
		else if (e[i] == 0)
		{
			for (int j = 0; j < m; j++)
			{
				aa[j].push_back(-1);
			}
		}

	}
	
	//处理d，添加松弛变量
	for (int i = 0; i < m; i++)
	{
		if (d[i] > 0)
		{
			for (int j = 0; j < m; j++)
			{
				if (j == i)
				{
					aa[j].push_back(-1);
				}
				else
					aa[j].push_back(0);
			}
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

void SimplexMatrix::readFile(const char* path)
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
	prepare();
}

void SimplexMatrix::print()
{
	cout << "m = " << m << endl;
	cout << "n = " << n << endl;
	cout << "c = ";
	for (int i = 0; i < n; i++)
	{
		cout << c[i] << "\t";
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

SimplexMatrix SimplexMatrix::transform()
{
	//另一个问题
	SimplexMatrix M;
	M.n = this->m + this->n;
	M.m = this->m;
	//b
	M.b = new double[M.m];
	for (int i = 0; i < M.m; i++)
	{
		M.b[i] = this->b[i];
	}
	//c
	M.c = new double[M.n];
	for (int i = 0; i < M.n; i++)
	{
		if (i < this->n)
			M.c[i] = 0;
		else
			M.c[i] = 1;
	}
	//a
	M.a = new double*[M.m];
	for (int i = 0; i < M.m; i++)
	{
		M.a[i] = new double[M.n];
		for (int j = 0; j < M.n; j++)
		{
			if (j < this->n)
				M.a[i][j] = this->a[i][j];
			else if (i == j - this->n)
				M.a[i][j] = 1;
			else
				M.a[i][j] = 0;
		}
	}
	//initial x
	M.x = new double[M.n];
	for (int i = 0; i < M.n; i++)
	{
		M.x[i] = 0;
	}
	//改变c
	for (int i = this->n; i < M.n; i++)
	{
		for (int j = 0; j < M.n; j++)
		{
			M.c[j] -= M.a[i - this->n][j];
		}
		M.z -= b[i - this->n];
	}
	return M;
}

void SimplexMatrix::print_root()
{
	if (x == NULL)
	{
		cout << "No root"<<endl;
		return;
	}
	//cout << "The root is: ";
	for (int i = 0; i < this->n; i++)
	{
		cout << x[i] << " ";
	}
	cout << endl;
}

int SimplexMatrix::do_counting()
{
	double **abc;  // 将成员变量a， b， c整合成单纯形表，以后的计算都应用于该表中
	double **n_abc;  // 更新单纯形表使用的临时矩阵
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

		in_base = -1;
		cur_min = 0;
		for (int i = 0; i < n; i++) {
			if (abc[m][i] < 0) {
				in_base = i;  // 寻找下表最小的负检验数为入基变量
				break;
			}
		}
		if (in_base < 0) {  // 找到最优解
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
			argmin_outbase = -1;
			min_outbase = DBL_MAX;
			for (int i = 0; i < m; i++) {  // 找出基变量
				if (abc[i][in_base] > 0) {
					temp = abc[i][n] / abc[i][in_base];
					//cout << i << '\t' << temp << '\t' << min_outbase << endl;
					if (temp < min_outbase) {
						min_outbase = temp;
						argmin_outbase = i;
					}
				}
			}

			if (argmin_outbase >= 0) {
				out_base = argmin_outbase;
			}
			else {
				//cout << "inbase" << in_base << endl;
				return 0;
			}
			for (int i = 0; i < m + 1; i++) {  // 转轴操作，更新单纯形表
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


SimplexMatrix SimplexMatrix::transport_back_ab()
{
	SimplexMatrix m;
	//m & n
	m.m = this->m;
	m.n = this->n - this->m;
	//b
	m.b = new double[m.m];
	for (int i = 0; i < m.m; i++)
	{
		m.b[i] = this->b[i];
	}
	//c
	m.c = new double[m.n];
	for (int i = 0; i < m.n; i++)
	{
		m.c[i] = 0;
	}
	//x
	m.x = new double[m.n];
	for (int i = 0; i < m.n; i++)
	{
		m.x[i] = 0;
	}
	//a
	m.a = new double*[m.m];
	for (int i = 0; i < m.m; i++)
	{
		m.a[i] = new double[m.n];
		for (int j = 0; j < m.n; j++)
		{
			m.a[i][j] = a[i][j];
		}
	}
	//z
	m.z = 0;
	return m;
}

void SimplexMatrix::copy_c(SimplexMatrix ori0)
{
	//拷贝原问题的c
	for (int i = 0; i < ori0.n; i++)
	{
		this->c[i] = ori0.c[i];
	}
}

void SimplexMatrix::make_normal(SimplexMatrix ori1) 
{
	int index;
	double factor;
	for (int i = 0; i < this->n; i++)
	{
		if (fabs(ori1.x[i]) > ZERO)
		{
			factor = c[i];
			for (int j = 0; j < this->m; j++)
			{
				if (this->a[j][i] == 1)
				{
					for (int k = 0; k < this->n; k++)
					{
						this->c[k] -= this->a[j][k] * factor;
					}
					this->z += factor * this->b[j];
					break;
				}
			}
			
		}
	}
}