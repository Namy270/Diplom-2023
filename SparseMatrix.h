#pragma once
#include <iostream>
#include <vector>
struct Col_Val
{
	int col;
	double val;
	Col_Val(int col, double val)
	{
		this->col = col;
		this->val = val;
	}
};

class SparseMatrix
{
public:
	std::vector<std::vector<Col_Val>> Matrix;
	int sizeV, sizeR;
	std::vector<int> csrR, csrC;
	std::vector<double> csrV;
	SparseMatrix(int sizeR);
	SparseMatrix();
	void add(int row, int col, double value);
	void resize(int size);
	void clear();
	void CSR();
	std::vector<double>& getVal();
	std::vector<int>& getColInd();
	std::vector<int>& getRowPtr();
	~SparseMatrix();
	//int size();
	//double& operator() (int row, int col);
};





