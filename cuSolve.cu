#include "cuSolve.cuh"
#include "device_launch_parameters.h"
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust\transform.h>
#include <thrust\functional.h>
#include <cuda_profiler_api.h>
#include <thrust/extrema.h>
#define Threads 32

__global__ void Poisk_Max(double* next_x, double* last_x, double* devMaxd_x)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    devMaxd_x[i] = abs(last_x[i] - next_x[i]);
}
double getNorm(double* dX, int n)
{
    /*thrust::device_ptr<double> d_ptrdX = thrust::device_pointer_cast(dX);
    double rez = *(thrust::max_element(d_ptrdX, d_ptrdX + n));*/
    double max = dX[0];
    for (int i = 0; i < n; i++)
    {
        if (dX[i] > max)
            max = dX[i];
    }

    return max;
}
__global__ void JacobiOnDevice(double* x_next, double* Value, int* Row, int* Col, double* x_now, double* b)
{
    // Optimization step 1: tiling
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j, start, N;
    double sigma;
    double diag;


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
    x_next[i] = (b[i] - sigma) / diag;
}
double* Solve_SLAY_device(std::vector<double>& Val, std::vector<int>& RowPtr, std::vector<int>& ColInd, std::vector<double>& vector_B, int iters)
{
    int n = RowPtr.size() - 1;
    double* result = new double[n];
    const double* csrValA = Val.data();
    const int* csrRowPtrA = RowPtr.data();
    const int* csrColIndA = ColInd.data();
    const double* B = vector_B.data();
    int RowSize = RowPtr.size();
    int ValSize = Val.size();

    double* Value_d, * b_d, * x_now_d, * x_next_d, * devMax_d_x;
    int* RowPtr_d, * ColInd_d;

    double* x_now, * x_next, * devMax;
    x_next = (double*)malloc(n * sizeof(double));
    x_now = (double*)malloc(n * sizeof(double));
    devMax = (double*)malloc(n * sizeof(double));

    for (int i = 0; i < n; i++)
        x_now[i] = x_next[i] = 0.0;

    // Allocate memory on the device
    cudaMalloc((void**)&Value_d, ValSize * sizeof(double));
    cudaMalloc((void**)&ColInd_d, ValSize * sizeof(int));
    cudaMalloc((void**)&RowPtr_d, RowSize * sizeof(int));
    cudaMalloc((void**)&b_d, n * sizeof(double));

    cudaMalloc((void**)&x_now_d, n * sizeof(double));
    cudaMalloc((void**)&x_next_d, n * sizeof(double));
    cudaMalloc((void**)&devMax_d_x, n * sizeof(double));


    // Copy data -> device
    cudaMemcpy(Value_d, csrValA, sizeof(double) * ValSize, cudaMemcpyHostToDevice);
    cudaMemcpy(ColInd_d, csrColIndA, sizeof(int) * ValSize, cudaMemcpyHostToDevice);
    cudaMemcpy(RowPtr_d, csrRowPtrA, sizeof(int) * RowSize, cudaMemcpyHostToDevice);
    cudaMemcpy(b_d, B, sizeof(double) * n, cudaMemcpyHostToDevice);
    cudaMemcpy(x_next_d, x_next, sizeof(double) * n, cudaMemcpyHostToDevice);
    cudaMemcpy(x_now_d, x_now, sizeof(double) * n, cudaMemcpyHostToDevice);

    dim3 threads = dim3(Threads);
    int Block = (int)ceil((double)n / Threads);
    dim3 blocks = dim3(Block);

    double dT1 = 1, Eps = 1e-8;
    int k = 0;

    //JacobiOnDevice << < blocks, threads >> > (x_now_d, Value_d, RowPtr_d, ColInd_d, x_next_d, b_d);
    //JacobiOnDevice << < blocks, threads >> > (x_next_d, Value_d, RowPtr_d, ColInd_d, x_now_d, b_d);
    //Poisk_Max << <blocks, threads >> > (x_next_d, x_now_d, devMax_d_x);
    //dT1 = getNorm(devMax_d_x, n);
    //JacobiOnDevice << < blocks, threads >> > (x_now_d, Value_d, RowPtr_d, ColInd_d, x_next_d, b_d);
    //Poisk_Max << <blocks, threads >> > (x_now_d, x_next_d, devMax_d_x);
    //dT2 = getNorm(devMax_d_x, n);
    //q = dT2 / dT1;
    //N = log(1e0 / Eps) / log(1e0 / q) + 1;

    for (k = 0;  k < iters; k++)
    {
        if (k % 2)
        {
            JacobiOnDevice << < blocks, threads >> > (x_now_d, Value_d, RowPtr_d, ColInd_d, x_next_d, b_d);
        }
        else
        {
            JacobiOnDevice << < blocks, threads >> > (x_next_d, Value_d, RowPtr_d, ColInd_d, x_now_d, b_d);
        }
       /* Poisk_Max << <blocks, threads >> > (x_next_d, x_now_d, devMax_d_x);
        cudaMemcpy(devMax, devMax_d_x, sizeof(double) * n, cudaMemcpyDeviceToHost);
        dT1 = getNorm(devMax, n);*/
    }
    cudaMemcpy(x_next, x_next_d, sizeof(double) * n, cudaMemcpyDeviceToHost);

    // Записываем результаты 
    for (int i = 0; i < n; i++)
        result[i] = x_next[i];

    // Освобождаем память 
    free(x_now);
    free(x_next);
    free(devMax);

    cudaFree(Value_d);
    cudaFree(ColInd_d);
    cudaFree(RowPtr_d);
    cudaFree(b_d);
    cudaFree(x_next_d);
    cudaFree(x_now_d);
    //std::cout << k << std::endl;
    return result;
}
