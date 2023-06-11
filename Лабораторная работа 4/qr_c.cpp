#include "mex.h"
#include <cmath>
#include <iostream>
#include <vector>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    mwSize m = mxGetM(prhs[0]);
    mwSize n = mxGetN(prhs[0]);
    
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(n, n, mxREAL);
    
    double* Q = mxGetPr(plhs[0]);
    double* R = mxGetPr(plhs[1]);
    const double* A = mxGetPr(prhs[0]);
    
    std::vector<double> v;
    v.resize(m);
    double sum;
    
    for (int j = 0; j < n; ++j){
        for (int x = 0; x < m; ++x){
            v[x] = A[x+j*m];
        }
        for (int i = 0; i <= j - 1; ++i){
            sum = 0;
            for (int k = 0; k < m; ++k){
                sum = sum + Q[k+i*m] * A[k+j*m];
            }
            for (int x = 0; x < m; ++x){
                v[x] = v[x] - sum * Q[x+i*m];
            }
        }
        sum = 0;
        for (int x = 0; x < m; ++x){
            sum += v[x] * v[x];
        }
        for (int x = 0; x < m; ++x){
            Q[x+j*m] = v[x] / std::sqrt(sum);
        } 
    }
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            R[i+j*n] = 0;
            for (int k = 0; k < m; k++)
                R[i+j*n] += Q[k+i*m] * A[k+j*m];
        }
    }
}