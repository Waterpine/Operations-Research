#pragma once
#include<iostream>
#include<string>
#include<vector>
#include<fstream>

#define MAX_M 500
#define ZERO 1e-4
using namespace std;

class SimplexMatrix
{
public:
	//构造函数
	SimplexMatrix() :n(0), m(0), z(0), e(NULL), d(NULL), b(NULL), c(NULL), a(NULL), x(NULL){};
	//读取问题，并且调用prepare进行预处理
	void readFile(const char* path);
	//打印a，b，c
	void print();
	//打印x
	void print_root();
	//由原问题得到二阶段法的问题
	SimplexMatrix transform();
	//返回二阶段法的最终表格的a和b
	SimplexMatrix transport_back_ab();
	//拷贝原问题的c
	void copy_c(SimplexMatrix ori0);
	//将拷贝的c由a，b和二阶段的根改为单纯形法能解的矩阵
	void make_normal(SimplexMatrix ori1);
	//单纯形法
	int do_counting();
	double z;
private:
	int n, m;
	int *e;
	int *d;
	double *b;
	double *c;
	double **a;
	double *x;
	//对读取的问题进行预处理
	void prepare();
};