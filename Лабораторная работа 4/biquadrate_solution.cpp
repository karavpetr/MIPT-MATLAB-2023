#include "mex.h"
#include <complex>
#include <iostream>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 3){
        mexErrMsgTxt("Three input arguments are required.");
    }
    int m1 = mxGetM(prhs[0]);
    int n1 = mxGetN(prhs[0]);
    int m2 = mxGetM(prhs[1]);
    int n2 = mxGetN(prhs[1]);
    int m3 = mxGetM(prhs[2]);
    int n3 = mxGetN(prhs[2]);
    if ((n1 != n2) || (n2 != n3) || (m1 != m2) || (m2 != m3)){
        mexErrMsgTxt("Input matrices are of different sizes.");
    }
    if (nlhs > 4){
        mexErrMsgTxt("Too many output arguments.");
    }

    int m = m1;
    int n = n1;

    plhs[0] = mxCreateDoubleMatrix(m, n, mxCOMPLEX);
    plhs[1] = mxCreateDoubleMatrix(m, n, mxCOMPLEX);
    plhs[2] = mxCreateDoubleMatrix(m, n, mxCOMPLEX);
    plhs[3] = mxCreateDoubleMatrix(m, n, mxCOMPLEX);

    double *x1r = mxGetPr(plhs[0]);
    double *x1i = mxGetPi(plhs[0]);
    double *x2r = mxGetPr(plhs[1]);
    double *x2i = mxGetPi(plhs[1]);
    double *x3r = mxGetPr(plhs[2]);
    double *x3i = mxGetPi(plhs[2]);
    double *x4r = mxGetPr(plhs[3]);
    double *x4i = mxGetPi(plhs[3]);

    double *ptrAr, *ptrAi, *ptrBr, *ptrBi, *ptrCr, *ptrCi;
    ptrAr = mxGetPr(prhs[0]);
    ptrAi = mxGetPi(prhs[0]);
    ptrBr = mxGetPr(prhs[1]);
    ptrBi = mxGetPi(prhs[1]);
    ptrCr = mxGetPr(prhs[2]);
    ptrCi = mxGetPi(prhs[2]);

    for (int i = 0; i < m*n; i++) {
        std::complex<double> a(ptrAr[i], ptrAi[i]);
        std::complex<double> b(ptrBr[i], ptrBi[i]);
        std::complex<double> c(ptrCr[i], ptrCi[i]);

        std::complex<double> first_sqrt = std::sqrt(b*b - 4.0*a*c);
        std::complex<double> x1 = std::sqrt((-b + first_sqrt) / 2.0 / a);
        std::complex<double> x2 = - std::sqrt((-b + first_sqrt) / 2.0 / a);
        std::complex<double> x3 = std::sqrt((-b - first_sqrt) / 2.0 / a);
        std::complex<double> x4 = - std::sqrt((-b - first_sqrt) / 2.0 / a);

        x1r[i] = x1.real();
        x1i[i] = x1.imag();
        x2r[i] = x2.real();
        x2i[i] = x2.imag();
        x3r[i] = x3.real();
        x3i[i] = x3.imag();
        x4r[i] = x4.real();
        x4i[i] = x4.imag();
    }
}
