/*
********��ż�����η�********
�������ڼ�����ȫ����0���޵�ʽԼ��
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
