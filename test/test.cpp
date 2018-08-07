#include<stdio.h>
#include<stdlib.h>
#include<time.h>

#define COL 500 
#define ROW 200 

int main(void) 
{
	FILE *fp;
	double M[ROW][COL];
	double C[COL]; // ��Ӧ���߱�����ϵ�� 
	int E[COL]; // ���ھ��߱������ŵ��ж� 
	int D[ROW]; // ����Լ���������ŵ��ж� 
	double B[ROW]; // ��ӦԼ�����Ҷ�ϵ�� 
	int i, j;
	int a, b;
	fp = fopen("dataB.txt","w"); // ���ļ�д��������� 
	int row = ROW;
	int col = COL;
	fprintf(fp,"%d ",COL); // д����߱������� 
	fprintf(fp,"%d\n",ROW); // д��Լ������ 
	for(i = 0; i < ROW; i++)
	{
		for(j = 0; j < COL; j++)
		{
			M[i][j] = 0;
		}
	}
	// ������ɾ��߱�����ϵ�� 
	for(i = 0; i < COL; i++)
	{
		C[i] = rand()/double(RAND_MAX) * 100;
		if(i == COL - 1)
		{
			fprintf(fp, "%lf\n",C[i]);
		}
		else
		{
			fprintf(fp, "%lf ",C[i]);
		}
	}
	// ��ʼ������һ�����н⣬��֤���Թ滮����н� 
	for(i = 0; i < ROW - 1; i++)
	{
		M[i][i] = 1; // 1��ʾ�ߵ���� 
		M[i+1][i] = -1; // -1��ʾ�ߵ��յ� 
	} 
	// �������ʣ�µı� 
	for(i = ROW - 1; i < COL; i++)
	{
		a = rand() % ROW;
		b = rand() % ROW;
		while(a == b)
		{
			a = rand() % ROW;
			b = rand() % ROW;
		}
		M[a][i] = 1;
		M[b][i] = -1;
	}
	for(i = 0; i < ROW; i++)
	{
		D[i] = 0;
		if(i == 0)
		{
			B[i] = 1;
		}
		else if(i == ROW - 1)
		{
			B[i] = -1;
		}
		else
		{
			B[i] = 0;
		}
	}
	for(i = 0; i < ROW; i++)
	{
		for(j = 0; j<COL; j++)
		{
			fprintf(fp,"%lf ",M[i][j]);
		}
		fprintf(fp, "%lf ", B[i]);
		fprintf(fp, "%d\n", D[i]);
	}
	for(i = 0; i < COL; i++)
	{
		E[i] = 1;
		if(i == COL - 1)
		{
			fprintf(fp, "%d\n",E[i]);
		}
		else
		{
			fprintf(fp, "%d ",E[i]);
		}
	}
}
