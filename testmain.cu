#ifdef TEST_INTEGRALS
#include <iostream>
int main() {
    const int numPairs = 2;
    ContractedIntegralPair h_pairs[numPairs];

    // integral 0
    for (int i = 0; i < numPairs; i++) {

        h_pairs[i].bra1.nPrim = 1;
        h_pairs[i].bra1.prims[0] = {0.5, 1.0, {0.0, 0.0, 0.0}, 0, 0, 0};

        h_pairs[i].bra2.nPrim = 1;
        h_pairs[i].bra2.prims[0] = {0.5, 1.0, {0.0, 0.0, 0.0}, 0, 0, 0};

        h_pairs[i].ket1.nPrim = 1;
        h_pairs[i].ket1.prims[0] = {0.5, 1.0, {0.0, 0.0, 0.0}, 0, 0, 0};

        h_pairs[i].ket2.nPrim = 1;
        h_pairs[i].ket2.prims[0] = {0.5, 1.0, {0.0, 0.0, 0.0}, 0, 0, 0};
    }

    for (int i = 0; i < 4; i++) {

        ContractedGaussian *cg = nullptr;
        if (i == 0) { // bra1: p_x function (lx=1)
            cg = &h_pairs[1].bra1;
            cg->nPrim = 2;
            cg->prims[0] = {0.4, 0.8, {0.0, 0.0, 0.0}, 1, 0, 0};
            cg->prims[1] = {0.2, 0.6, {0.0, 0.0, 0.0}, 1, 0, 0};
        } else if (i == 1) { // bra2: s function (l=0)
            cg = &h_pairs[1].bra2;
            cg->nPrim = 2;
            cg->prims[0] = {0.4, 0.8, {0.1, 0.0, 0.0}, 0, 0, 0};
            cg->prims[1] = {0.2, 0.6, {0.1, 0.0, 0.0}, 0, 0, 0};
        } else if (i == 2) { // ket1: p_y function (ly=1)
            cg = &h_pairs[1].ket1;
            cg->nPrim = 2;
            cg->prims[0] = {0.4, 0.8, {0.0, 0.1, 0.0}, 0, 1, 0};
            cg->prims[1] = {0.2, 0.6, {0.0, 0.1, 0.0}, 0, 1, 0};
        } else if (i == 3) { 
            cg = &h_pairs[1].ket2;
            cg->nPrim = 2;
            cg->prims[0] = {0.4, 0.8, {0.0, 0.0, 0.1}, 0, 0, 0};
            cg->prims[1] = {0.2, 0.6, {0.0, 0.0, 0.1}, 0, 0, 0};
        }
    }

    double h_results[numPairs] = {0.0};

    // eval
    evaluate_contracted_two_electron_integrals(h_pairs, h_results, numPairs);

    for (int i = 0; i < numPairs; i++) {
        std::cout << "Contracted Integral " << i << ": " << h_results[i] << std::endl;
    }

    return 0;
}
#endif
