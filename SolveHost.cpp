#include "SolveHost.h"
#include  <omp.h>
double getNorm(double* x_last, double* x_next, int n)
{
    double max = -1;
    double razn;
    for (int i = 0; i < n; i++)
    {
        razn = abs(x_last[i] - x_next[i]);
        if (razn > max)
            max = razn;
    }
    return max;
}
void jacobiOnHost(std::vector<double>& x_next, const double* Value, const int* Row, const int* Col, std::vector<double>& x_now, const double* b, int Ni)
{
    int j, start, N;
    double sigma;
    double diag = 1;
    #pragma omp parallel for
    for (int i = 0; i < Ni; i++)
    {
        start = Row[i];
        N = Row[i + 1] - Row[i];
        sigma = 0.0;
        for (int l = 0; l < N; l++)
        {
            j = Col[start + l];
            if (i != j)
                sigma += x_now[j] * Value[start + l];
            else
                diag = Value[start + l];
        }
        double help = x_next[i];
        x_next[i] = (b[i] - sigma) / diag;
    }
}
double* Solve_SLAY_host(std::vector<double>& Val, std::vector<int>& RowPtr, std::vector<int>& ColInd, std::vector<double>& vector_B, int iters)
{
    int n = RowPtr.size() - 1;
    double* result = new double[n];
    const double* csrValA = Val.data();
    const int* csrRowPtrA = RowPtr.data();
    const int* csrColIndA = ColInd.data();
    const double* B = vector_B.data();
    int RowSize = RowPtr.size();
    int ValSize = Val.size();
    std::vector<double> x_next(n), x_now(n);
    for (int i = 0; i < n; i++)
        x_now[i] = x_next[i] = 0.0;
    double devMax_x;
    double dX = 1;
    int k = 0;
    //#pragma omp parallel for
    for (k = 0; k < iters; k++)
    {
        /* if (k % 2)
             jacobiOnHost(x_now, csrValA, csrRowPtrA, csrColIndA, x_next, B, n);
         else*/
        jacobiOnHost(x_next, csrValA, csrRowPtrA, csrColIndA, x_now, B, n);
        //dX = getNorm(x_next.data(), x_now.data(), n);
        x_now = x_next;
    }


    // Записываем результаты 
    for (int i = 0; i < n; i++)
        result[i] = x_next[i];
    //std::cout << k << std::endl;
    // Освобождаем память 
    return result;
}
