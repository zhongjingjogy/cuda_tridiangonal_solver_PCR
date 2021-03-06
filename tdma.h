#ifndef TDMA_H
#define TDMA_H
/*
Tri-Diagonal Matrix Algorithm for one dimension linear equation
link: https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
@ a: lower array of the tri-diagonal matrix
@ b: exact diagonal array of the tri-diagonal matrix
@ c: upper array of the tri-diagonal matrix
@ rhs: right hand side of the linear equations
@ old: solution to linear equations
@ Nx: dimension of the equations
*/
template <typename Type>
void TDMA(Type *a, Type *b, Type *c, Type *rhs, Type *old, int Nx) {
    Type *beta = new Type[Nx];
    Type *y = new Type[Nx];
    beta[0] = c[0] / b[0];
    for (int nx = 1; nx < Nx - 1; nx++) {
        beta[nx] = c[nx] / (b[nx] - a[nx] * beta[nx - 1]);
    }
    y[0] = rhs[0] / b[0];
    for (int nx = 1; nx < Nx; nx++) {
        y[nx] = (rhs[nx] - a[nx] * y[nx - 1]) / (b[nx] - a[nx] * beta[nx - 1]);
    }
    old[Nx - 1] = y[Nx - 1];
    for (int nx = Nx - 2; nx >= 0; nx--) {
        old[nx] = y[nx] - beta[nx] * old[nx + 1];
    }
    delete[] beta;
    delete[] y;
}

template <typename Type>
void TDMA(std::vector<Type> &a, std::vector<Type> &b, std::vector<Type> &c,
          std::vector<Type> &rhs, std::vector<Type> &old, int Nx) {
    Type *beta = new Type[Nx];
    Type *y = new Type[Nx];
    beta[0] = c[0] / b[0];
    for (int nx = 1; nx < Nx - 1; nx++) {
        beta[nx] = c[nx] / (b[nx] - a[nx] * beta[nx - 1]);
    }
    y[0] = rhs[0] / b[0];
    for (int nx = 1; nx < Nx; nx++) {
        y[nx] = (rhs[nx] - a[nx] * y[nx - 1]) / (b[nx] - a[nx] * beta[nx - 1]);
    }
    old[Nx - 1] = y[Nx - 1];
    for (int nx = Nx - 2; nx >= 0; nx--) {
        old[nx] = y[nx] - beta[nx] * old[nx + 1];
    }
    delete[] beta;
    delete[] y;
}

template <typename Type>
void TDMAsolve(Type *a, Type *b, Type *c, Type *d, int n) {
    /*
    // n is the number of unknowns

    |b0 c0 0 ||x0| |d0|
    |a1 b1 c1||x1|=|d1|
    |0  a2 b2||x2| |d2|

    1st iteration: b0x0 + c0x1 = d0 -> x0 + (c0/b0)x1 = d0/b0 ->

        x0 + g0x1 = r0               where g0 = c0/b0        , r0 = d0/b0

    2nd iteration:     | a1x0 + b1x1   + c1x2 = d1
        from 1st it.: -| a1x0 + a1g0x1        = a1r0
                    -----------------------------
                          (b1 - a1g0)x1 + c1x2 = d1 - a1r0

        x1 + g1x2 = r1               where g1=c1/(b1 - a1g0) , r1 = (d1 -
    a1r0)/(b1 - a1g0)

    3rd iteration:      | a2x1 + b2x2   = d2
        from 2nd it. : -| a2x1 + a2g1x2 = a2r2
                       -----------------------
                       (b2 - a2g1)x2 = d2 - a2r2
        x2 = r2                      where                     r2 = (d2 -
    a2r2)/(b2 - a2g1) Finally we have a triangular matrix: |1  g0 0 ||x0| |r0|
    |0  1  g1||x1|=|r1|
    |0  0  1 ||x2| |r2|

    Condition: ||bi|| > ||ai|| + ||ci||

    in this version the c matrix reused instead of g
    and             the d matrix reused instead of r and x matrices to report
    results Written by Keivan Moradi, 2014
    */
    n--;  // since we start from x0 (not x1)
    c[0] /= b[0];
    d[0] /= b[0];

    for (int i = 1; i < n; i++) {
        c[i] /= b[i] - a[i] * c[i - 1];
        d[i] = (d[i] - a[i] * d[i - 1]) / (b[i] - a[i] * c[i - 1]);
    }

    d[n] = (d[n] - a[n] * d[n - 1]) / (b[n] - a[n] * c[n - 1]);

    for (int i = n; i-- > 0;) {
        d[i] -= c[i] * d[i + 1];
    }
}

// from:
// https://en.wikibooks.org/wiki/Algorithm_Implementation/Linear_Algebra/Tridiagonal_matrix_algorithm
template <typename Type>
void TDMAsolve(std::vector<Type> &a, std::vector<Type> &b, std::vector<Type> &c,
               std::vector<Type> &d, int n) {
    /*
    // n is the number of unknowns

    |b0 c0 0 ||x0| |d0|
    |a1 b1 c1||x1|=|d1|
    |0  a2 b2||x2| |d2|

    1st iteration: b0x0 + c0x1 = d0 -> x0 + (c0/b0)x1 = d0/b0 ->

        x0 + g0x1 = r0               where g0 = c0/b0        , r0 = d0/b0

    2nd iteration:     | a1x0 + b1x1   + c1x2 = d1
        from 1st it.: -| a1x0 + a1g0x1        = a1r0
                    -----------------------------
                          (b1 - a1g0)x1 + c1x2 = d1 - a1r0

        x1 + g1x2 = r1               where g1=c1/(b1 - a1g0) , r1 = (d1 -
    a1r0)/(b1 - a1g0)

    3rd iteration:      | a2x1 + b2x2   = d2
        from 2nd it. : -| a2x1 + a2g1x2 = a2r2
                       -----------------------
                       (b2 - a2g1)x2 = d2 - a2r2
        x2 = r2                      where                     r2 = (d2 -
    a2r2)/(b2 - a2g1) Finally we have a triangular matrix: |1  g0 0 ||x0| |r0|
    |0  1  g1||x1|=|r1|
    |0  0  1 ||x2| |r2|

    Condition: ||bi|| > ||ai|| + ||ci||

    in this version the c matrix reused instead of g
    and             the d matrix reused instead of r and x matrices to report
    results Written by Keivan Moradi, 2014
    */
    n--;  // since we start from x0 (not x1)
    c[0] /= b[0];
    d[0] /= b[0];

    for (int i = 1; i < n; i++) {
        c[i] /= b[i] - a[i] * c[i - 1];
        d[i] = (d[i] - a[i] * d[i - 1]) / (b[i] - a[i] * c[i - 1]);
    }

    d[n] = (d[n] - a[n] * d[n - 1]) / (b[n] - a[n] * c[n - 1]);

    for (int i = n; i-- > 0;) {
        d[i] -= c[i] * d[i + 1];
    }
}

#endif