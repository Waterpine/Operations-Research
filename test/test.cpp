#include<stdio.h>
#include<stdlib.h>
#include<time.h>

#define COL 500 
#define ROW 200 

int main(void) 
{
	FILE *fp;
	double M[ROW][COL];
	double C[COL]; // 对应决策变量的系数 
	int E[COL]; // 对于决策变量符号的判断 
	int D[ROW]; // 对于约束条件符号的判断 
	double B[ROW]; // 对应约束的右端系数 
	int i, j;
	int a, b;
	fp = fopen("dataB.txt","w"); // 打开文件写入测试数据 
	int row = ROW;
	int col = COL;
	fprintf(fp,"%d ",COL); // 写入决策变量个数 
	fprintf(fp,"%d\n",ROW); // 写入约束个数 
	for(i = 0; i < ROW; i++)
	{
		for(j = 0; j < COL; j++)
		{
			M[i][j] = 0;
		}
	}
	// 随机生成决策变量的系数 
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
	// 初始化生成一个可行解，保证线性规划务必有解 
	for(i = 0; i < ROW - 1; i++)
	{
		M[i][i] = 1; // 1表示边的起点 
		M[i+1][i] = -1; // -1表示边的终点 
	} 
	// 随机生成剩下的边 
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
