#include <cmath>
#include <cstdio>
#include <cuda_runtime.h>
#include <data.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846


__device__ double dist2(double3 a, double3 b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    double dz = a.z - b.z;
    return dx*dx + dy*dy + dz*dz;
}

// F_m(T) = ∫₀¹ u^(2m) exp(-T u^2) du
__device__ double boys_function(int m, double T) {
    if (T < 1e-8) {
        return 1.0/(2*m+1);
    } else {
        double F0 = 0.5 * sqrt(M_PI/T) * erf(sqrt(T));
        double F = F0;
        for (int i = 1; i <= m; i++) {
            F = ((2.0*i - 1)*F - exp(-T))/(2.0*T);
        }
        return F;
    }
}

__device__ void computeHermiteCoeffs(int la, int lb, double P, double base,
                                     double A, double B, double p, double *E_out) {
    int size_i = la + 1;
    int size_j = lb + 1;
    int total_size = size_i * size_j;
    // Temporary array F[i][j] stored in 1D (assume total_size <= 100)
    double F[100];
    for (int i = 0; i < total_size; i++) {
        F[i] = 0.0;
    }
    //base case
    F[0] = base;
    
    for (int i = 1; i < size_i; i++) {
        double term1 = (P - A) * F[(i-1)*size_j + 0];
        double term2 = 0.0;
        if (i - 1 >= 1) {
            term2 = ((double)(i - 1))/(2.0 * p) * F[(i-2)*size_j + 0];
        }
        F[i*size_j + 0] = term1 + term2;
    }
    
    for (int j = 1; j < size_j; j++) {
        double term1 = (P - B) * F[0*size_j + (j-1)];
        double term2 = 0.0;
        if (j - 1 >= 1) {
            term2 = ((double)(j - 1))/(2.0 * p) * F[0*size_j + (j-2)];
        }
        F[0*size_j + j] = term1 + term2;
    }
    // recurrence
    for (int i = 1; i < size_i; i++) {
        for (int j = 1; j < size_j; j++) {
            double term1 = (P - A) * F[(i-1)*size_j + j];
            double term2 = (P - B) * F[i*size_j + (j-1)];
            double term3 = ((double)(i+j-1))/(2.0 * p) * F[(i-1)*size_j + (j-1)];
            F[i*size_j + j] = term1 + term2 + term3;
        }
    }
    // sum reccurence
    int L = la + lb;
    for (int t = 0; t <= L; t++) {
        double sum = 0.0;
        for (int i = 0; i <= la; i++) {
            int j = t - i;
            if (j < 0 || j > lb) continue;
            sum += F[i*size_j + j];
        }
        E_out[t] = sum;
    }
}
