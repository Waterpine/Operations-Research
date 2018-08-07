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
	//���캯��
	SimplexMatrix() :n(0), m(0), z(0), e(NULL), d(NULL), b(NULL), c(NULL), a(NULL), x(NULL){};
	//��ȡ���⣬���ҵ���prepare����Ԥ����
	void readFile(const char* path);
	//��ӡa��b��c
	void print();
	//��ӡx
	void print_root();
	//��ԭ����õ����׶η�������
	SimplexMatrix transform();
	//���ض��׶η������ձ���a��b
	SimplexMatrix transport_back_ab();
	//����ԭ�����c
	void copy_c(SimplexMatrix ori0);
	//��������c��a��b�Ͷ��׶εĸ���Ϊ�����η��ܽ�ľ���
	void make_normal(SimplexMatrix ori1);
	//�����η�
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
	//�Զ�ȡ���������Ԥ����
	void prepare();
};