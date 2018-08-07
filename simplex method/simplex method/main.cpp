/*
****************单纯形法*****************
修改readFile中的参数即可读取文件并给出解
约束上限为500
*/
#include"simplex.h"
#include<cmath>
#pragma warning(disable:4996)
int main()
{
	SimplexMatrix M;
	//读取并且预处理
	M.readFile("临界.txt");
	SimplexMatrix M1;
	//转化为两步法的问题
	M1 = M.transform();
	//第一次迭代
	M1.do_counting();
	//无解
	if (M1.z > ZERO)
	{
		cout << -1 << endl;
	}
	//有解
	else
	{
		SimplexMatrix M2;
		//将矩阵返回
		M2 = M1.transport_back_ab();
		M2.copy_c(M);
		M2.make_normal(M1);
		//无限最优解
		if (M2.do_counting() == 0)
		{
			cout << 0 << endl;
			return 0;
		}
		//有限最优解
		cout << 1 << endl;
		cout << M2.z << endl;
		M2.print_root();
	}
}
