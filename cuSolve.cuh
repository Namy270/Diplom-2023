#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <cuda_runtime.h>
#include <cusparse.h>
#include <vector>
double* Solve_SLAY_device(std::vector<double>& Val, std::vector<int>& RowPtr, std::vector<int>& ColInd, std::vector<double>& vector_B, int iters);
