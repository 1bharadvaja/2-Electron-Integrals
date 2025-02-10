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

//will optimize these recurrences later on


__device__ double compute_eri_angular(const GaussianPrimitive &a, const GaussianPrimitive &b,
                                        const GaussianPrimitive &c, const GaussianPrimitive &d) {
    // bra pair (ab)
    double alpha = a.exponent;
    double beta  = b.exponent;
    double p = alpha + beta;
    double3 A = a.center;
    double3 B = b.center;
    double3 P;
    P.x = (alpha * A.x + beta * B.x) / p;
    P.y = (alpha * A.y + beta * B.y) / p;
    P.z = (alpha * A.z + beta * B.z) / p;
    
    // ket pair (cd)
    double gamma = c.exponent;
    double delta = d.exponent;
    double q = gamma + delta;
    double3 C = c.center;
    double3 D = d.center;
    double3 Q;
    Q.x = (gamma * C.x + delta * D.x) / q;
    Q.y = (gamma * C.y + delta * D.y) / q;
    Q.z = (gamma * C.z + delta * D.z) / q;

    double base_bra_x = exp( - (alpha * beta / p) * ((A.x - B.x)*(A.x - B.x)) );
    double base_bra_y = exp( - (alpha * beta / p) * ((A.y - B.y)*(A.y - B.y)) );
    double base_bra_z = exp( - (alpha * beta / p) * ((A.z - B.z)*(A.z - B.z)) );
    
    double base_ket_x = exp( - (gamma * delta / q) * ((C.x - D.x)*(C.x - D.x)) );
    double base_ket_y = exp( - (gamma * delta / q) * ((C.y - D.y)*(C.y - D.y)) );
    double base_ket_z = exp( - (gamma * delta / q) * ((C.z - D.z)*(C.z - D.z)) );
    
    // angular momenta
    // Bra pair: orbital a and b.
    int la_x = a.lx, la_y = a.ly, la_z = a.lz;
    int lb_x = b.lx, lb_y = b.ly, lb_z = b.lz;

    // Ket pair: orbital c and d.
    int lc_x = c.lx, lc_y = c.ly, lc_z = c.lz;
    int ld_x = d.lx, ld_y = d.ly, ld_z = d.lz;
    
    // hermite coefficients for bra pair
    int nbra_x = la_x + lb_x;  // maximum order in x is la_x+lb_x
    int nbra_y = la_y + lb_y;
    int nbra_z = la_z + lb_z;
    const int maxOrder = 10; // (assumes angular momenta are small)
    double Ebra_x[maxOrder] = {0.0};
    double Ebra_y[maxOrder] = {0.0};
    double Ebra_z[maxOrder] = {0.0};
    computeHermiteCoeffs(la_x, lb_x, P.x, base_bra_x, A.x, B.x, p, Ebra_x);
    computeHermiteCoeffs(la_y, lb_y, P.y, base_bra_y, A.y, B.y, p, Ebra_y);
    computeHermiteCoeffs(la_z, lb_z, P.z, base_bra_z, A.z, B.z, p, Ebra_z);
    
    // hermite coefficients for ket 
    int nket_x = lc_x + ld_x;
    int nket_y = lc_y + ld_y;
    int nket_z = lc_z + ld_z;
    double Eket_x[maxOrder] = {0.0};
    double Eket_y[maxOrder] = {0.0};
    double Eket_z[maxOrder] = {0.0};
    computeHermiteCoeffs(lc_x, ld_x, Q.x, base_ket_x, C.x, D.x, q, Eket_x);
    computeHermiteCoeffs(lc_y, ld_y, Q.y, base_ket_y, C.y, D.y, q, Eket_y);
    computeHermiteCoeffs(lc_z, ld_z, Q.z, base_ket_z, C.z, D.z, q, Eket_z);
    
    // combine pairs
    int Lx = nbra_x + nket_x;
    int Ly = nbra_y + nket_y;
    int Lz = nbra_z + nket_z;
    double Fx[20] = {0.0}; // allocate arrays large enough (here assume L < 20)
    double Fy[20] = {0.0};
    double Fz[20] = {0.0};
    for (int i = 0; i <= nbra_x; i++) {
        for (int j = 0; j <= nket_x; j++) {
            Fx[i+j] += Ebra_x[i] * Eket_x[j];
        }
    }
    for (int i = 0; i <= nbra_y; i++) {
        for (int j = 0; j <= nket_y; j++) {
            Fy[i+j] += Ebra_y[i] * Eket_y[j];
        }
    }
    for (int i = 0; i <= nbra_z; i++) {
        for (int j = 0; j <= nket_z; j++) {
            Fz[i+j] += Ebra_z[i] * Eket_z[j];
        }
    }
    
    // --- Compute Boys function argument --- 
    double dx = P.x - Q.x;
    double dy = P.y - Q.y;
    double dz = P.z - Q.z;
    double T = (p * q / (p + q)) * (dx*dx + dy*dy + dz*dz);
    
    // --- Triple sum over Cartesian indices ---
    double sum = 0.0;
    for (int tx = 0; tx <= Lx; tx++) {
        for (int ty = 0; ty <= Ly; ty++) {
            for (int tz = 0; tz <= Lz; tz++) {
                int order = tx + ty + tz;
                double boys = boys_function(order, T);
                sum += Fx[tx] * Fy[ty] * Fz[tz] * boys;
            }
        }
    }
    
    double prefactor = 2.0 * pow(M_PI, 2.5) / (p * q * sqrt(p + q));
    double result = a.coefficient * b.coefficient * c.coefficient * d.coefficient *
                    prefactor * sum;
    return result;
}


__device__ double compute_contracted_eri(const ContractedGaussian &A,
                                           const ContractedGaussian &B,
                                           const ContractedGaussian &C,
                                           const ContractedGaussian &D) {
    double sum = 0.0;
    for (int i = 0; i < A.nPrim; i++) {
        for (int j = 0; j < B.nPrim; j++) {
            for (int k = 0; k < C.nPrim; k++) {
                for (int l = 0; l < D.nPrim; l++) {

                    double term = compute_eri_angular(A.prims[i],
                                                      B.prims[j],
                                                      C.prims[k],
                                                      D.prims[l]);
                    sum += term;
                }
            }
        }
    }
    return sum;
}

// kernel
__global__ void compute_contracted_eri_kernel(const ContractedIntegralPair* pairs,
                                                double* results,
                                                int numPairs) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < numPairs) {
        const ContractedIntegralPair &pair = pairs[idx];
        double eri = compute_contracted_eri(pair.bra1, pair.bra2,
                                            pair.ket1, pair.ket2);
        results[idx] = eri;
    }
}

//API Functions!
extern "C" void evaluate_contracted_two_electron_integrals(const ContractedIntegralPair* h_pairs,
                                                             double* h_results,
                                                             int numPairs) {
    ContractedIntegralPair* d_pairs = nullptr;
    double* d_results = nullptr;
    size_t pairsSize   = numPairs * sizeof(ContractedIntegralPair);
    size_t resultsSize = numPairs * sizeof(double);

    cudaError_t err;
    err = cudaMalloc((void**)&d_pairs, pairsSize);
    if (err != cudaSuccess) {
        fprintf(stderr, "cudaMalloc for d_pairs failed: %s\n", cudaGetErrorString(err));
        return;
    }
    err = cudaMalloc((void**)&d_results, resultsSize);
    if (err != cudaSuccess) {
        fprintf(stderr, "cudaMalloc for d_results failed: %s\n", cudaGetErrorString(err));
        cudaFree(d_pairs);
        return;
    }
    err = cudaMemcpy(d_pairs, h_pairs, pairsSize, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy to device failed: %s\n", cudaGetErrorString(err));
        cudaFree(d_pairs);
        cudaFree(d_results);
        return;
    }
    int threadsPerBlock = 256;
    int blocks = (numPairs + threadsPerBlock - 1) / threadsPerBlock;
    compute_contracted_eri_kernel<<<blocks, threadsPerBlock>>>(d_pairs, d_results, numPairs);
    err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "Kernel launch error: %s\n", cudaGetErrorString(err));
        cudaFree(d_pairs);
        cudaFree(d_results);
        return;
    }
    err = cudaMemcpy(h_results, d_results, resultsSize, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy from device failed: %s\n", cudaGetErrorString(err));
    }
    cudaFree(d_pairs);
    cudaFree(d_results);
}


