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

