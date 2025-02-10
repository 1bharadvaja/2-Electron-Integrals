#include <cmath>
#include <cstdio>
#include <cuda_runtime.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ---------------------------------------------------------------------------
// Data Structures
// ---------------------------------------------------------------------------

// defines Gaussian Primitive by its exponent, coeff, center, and ang. momentum in 3 directions
struct GaussianPrimitive {
    double exponent;      // Gaussian exponent α
    double coefficient;   // Contraction coefficient (including normalization)
    double3 center;       // Center of the Gaussian (x,y,z)
    int lx, ly, lz;       // Angular momentum in x, y, and z directions
};

// A quartet of primitives 
struct IntegralPair {
    GaussianPrimitive bra1;  // orbital a (first in bra pair)
    GaussianPrimitive bra2;  // orbital b (second in bra pair)
    GaussianPrimitive ket1;  // orbital c (first in ket pair)
    GaussianPrimitive ket2;  // orbital d (second in ket pair)
};
