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
void TDMA(float *a, float *b, float *c, float *rhs, float *old, int Nx) {
    float *beta = new float[Nx];
    float *y = new float[Nx];
    beta[0] = c[0] / b[0];
    for(int nx=1; nx<Nx-1; nx++) {
        beta[nx] = c[nx] / (b[nx] - a[nx]*beta[nx-1]);
    }
    y[0] = rhs[0] / b[0];
    for(int nx=1; nx<Nx; nx++) {
        y[nx] = (rhs[nx] - a[nx]*y[nx-1]) / (b[nx] - a[nx]*beta[nx-1]);
    }
    old[Nx-1] = y[Nx-1];
    for(int nx=Nx-2; nx>=0; nx--) {
        old[nx] = y[nx] - beta[nx]*old[nx+1];
    }
    delete []beta;
    delete []y;
}

#endif