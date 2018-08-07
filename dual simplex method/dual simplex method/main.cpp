/*
********对偶单纯形法********
仅适用于检验数全大于0且无等式约束
*/
#include "dualsimplex.h"
#include<cmath>
#pragma warning(disable:4996)
#define ZERO 1e-10

int main()
{
	char ch;
	DualSimplexMatrix M;
	M.readFile("data.txt");
	M.prepare();
	//M.print();
	//cout << "****************************************************\n";
	if (M.do_counting() == 1)
	{
		cout << 1 << endl;
		cout << M.z << endl;
		M.print_root();
		
	}
	else
	{
		cout << 0 << endl;
	}
	scanf("%c", &ch);
}
