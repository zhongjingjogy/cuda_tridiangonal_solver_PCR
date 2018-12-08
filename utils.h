#ifndef UTILS_H
#define UTILS_H
double maxarr(double *arr, int Nx) {
    double v = arr[0];
    for (int nx = 0; nx < Nx; nx++) {
        if (v < arr[nx]) {
            v = arr[nx];
        }
    }
    return v;
}

double minarr(double *arr, int Nx) {
    double v = arr[0];
    for (int nx = 0; nx < Nx; nx++) {
        if (v > arr[nx]) {
            v = arr[nx];
        }
    }
    return v;
}

double averarr(double *arr, int Nx) {
    double v = 0.0;
    for (int nx = 0; nx < Nx; nx++) {
        v += arr[nx];
    }
    return v / double(Nx);
}

template <typename Type>
Type CpuResultCheck(Type *a, Type *b, Type *c, Type *rhs, Type *old, int Nx) {
    Type v = -1.0;
    Type v1 = fabs(b[0] * old[0] + c[0] * old[1] - rhs[0]);
    if (v < v1) {
        v = v1;
    }
    for (int nx = 1; nx < Nx - 1; nx++) {
        v1 = fabs(a[nx] * old[nx - 1] + b[nx] * old[nx] + c[nx] * old[nx + 1] -
                  rhs[nx]);
        if (v < v1) {
            v = v1;
        }
    }
    v1 = fabs(b[Nx - 1] * old[Nx - 1] + a[Nx - 1] * old[Nx - 2] - rhs[Nx - 1]);
    if (v < v1) {
        v = v1;
    }

    return v;
}

template <typename Type>
Type CpuResultCheck(std::vector<Type> &a, std::vector<Type> &b,
                      std::vector<Type> &c, std::vector<Type> &rhs,
                      std::vector<Type> &old, int Nx) {
    Type v = -1.0;
    Type v1 = fabs(b[0] * old[0] + c[0] * old[1] - rhs[0]);
    if (v < v1) {
        v = v1;
    }
    for (int nx = 1; nx < Nx - 1; nx++) {
        v1 = fabs(a[nx] * old[nx - 1] + b[nx] * old[nx] + c[nx] * old[nx + 1] -
                  rhs[nx]);
        if (v < v1) {
            v = v1;
        }
    }
    v1 = fabs(b[Nx - 1] * old[Nx - 1] + a[Nx - 1] * old[Nx - 2] - rhs[Nx - 1]);
    if (v < v1) {
        v = v1;
    }

    return v;
}

int correctcount(bool *result, int trynumber) {
    int count = 0;
    for (int n = 0; n < trynumber; n++) {
        if (result[n]) {
            count++;
        }
    }
    return count;
}
#endif