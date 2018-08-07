#pragma once
#include<iostream>
#include<string>
#include<vector>
#include<fstream>

#define MAX_LINE 2048
#define MAX_M 2000
using namespace std;

class DualSimplexMatrix
{
public:
	void readFile(const char* path);
	void print();
	void print_root(); 
	int do_counting();
	//对读取的问题进行预处理
	void prepare();
	void dual();
	void make_normal();
	double z;
private:
	int n, m;
	int *e;
	int *d;
	double *b;
	double *c;
	double **a;
	double **dual_a;
	double *x;
};