/*
****************�����η�*****************
�޸�readFile�еĲ������ɶ�ȡ�ļ���������
Լ������Ϊ500
*/
#include"simplex.h"
#include<cmath>
#pragma warning(disable:4996)
int main()
{
	SimplexMatrix M;
	//��ȡ����Ԥ����
	M.readFile("�ٽ�.txt");
	SimplexMatrix M1;
	//ת��Ϊ������������
	M1 = M.transform();
	//��һ�ε���
	M1.do_counting();
	//�޽�
	if (M1.z > ZERO)
	{
		cout << -1 << endl;
	}
	//�н�
	else
	{
		SimplexMatrix M2;
		//�����󷵻�
		M2 = M1.transport_back_ab();
		M2.copy_c(M);
		M2.make_normal(M1);
		//�������Ž�
		if (M2.do_counting() == 0)
		{
			cout << 0 << endl;
			return 0;
		}
		//�������Ž�
		cout << 1 << endl;
		cout << M2.z << endl;
		M2.print_root();
	}
}
