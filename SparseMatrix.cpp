#include "SparseMatrix.h"

SparseMatrix::SparseMatrix(int size)
{
	Matrix.resize(size);
	sizeV = 0;
	sizeR = size + 1;
}

SparseMatrix::SparseMatrix()
{
	sizeV = 0;
	sizeR = 0;
}

void SparseMatrix::add(int row, int col, double value)
{
	Matrix[row].push_back(Col_Val(col, value));
	sizeV++;
}

void SparseMatrix::resize(int size)
{
	SparseMatrix::clear();
	Matrix.resize(size);
	sizeR = size + 1;
	sizeV = 0;
}

void SparseMatrix::CSR()
{
	csrR.push_back(0);
	for (int i = 0; i < Matrix.size(); i++)
	{
		for (int j = 0; j < Matrix[i].size(); j++)
		{
			csrV.push_back(Matrix[i][j].val);
			csrC.push_back(Matrix[i][j].col);
		}
		csrR.push_back(csrC.size());
	}
}

std::vector<double>& SparseMatrix::getVal()
{
	return csrV;
}

std::vector<int>& SparseMatrix::getColInd()
{
	return csrC;
}

std::vector<int>& SparseMatrix::getRowPtr()
{
	return csrR;
}

void SparseMatrix::clear()
{
	for (int i = 0; i < Matrix.size(); i++)
		Matrix[i].clear();
	Matrix.clear();
	csrR.clear();
	csrC.clear();
	csrV.clear();
	sizeR = sizeV = 0;
}
SparseMatrix::~SparseMatrix()
{
	clear();
}
//int MatrixSparse::size()
//{
//	return sizeR - 1;
//}

//double& MatrixSparse::operator()(int row, int col)
//{
//	Row[row]++;
//	Col.push_back(col);
//	Val.push_back(0);
//	sizeV++;
//	return Val.back();
//}