#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "pcr/PCR_Class.h"
#include "tdma/tdma.h"

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

bool CpuResultCheck(float *a, float *b, float *c, float *rhs, float *old,
                    int Nx) {
    if (fabs(b[0] * old[0] + c[0] * old[1] - rhs[0]) > 1.0e-6) {
        return false;
    }
    for (int nx = 1; nx < Nx - 1; nx++) {
        if (fabs(a[nx] * old[nx - 1] + b[nx] * old[nx] +
                 c[n] * old[nx + 1] - rhs[nx]) > 1.0e-6) {
            return false;
        }
    }
    if (fabs(b[Nx - 1] * old[Nx - 1] + a[Nx - 1] * old[Nx - 2] - rhs[Nx - 1]) >
        1.0e-6) {
        return false;
    }
}

bool GpuResultCheck(float *alist, float *blist, float *clist, float *dlist,
                    float *xlist, int DMax) {
    if (fabs(blist[0] * xlist[0] + clist[0] * xlist[1] - dlist[0]) > 1.0e-6) {
        return false;
    }
    for (int nx = 1; nx < DMax - 1; nx++) {
        if (fabs(alist[nx] * xlist[nx - 1] + blist[nx] * xlist[nx] +
                 clist[nx] * xlist[nx + 1] - dlist[nx]) > 1.0e-6) {
            return false;
        }
    }
    if (fabs(blist[Nx - 1] * xlist[Nx - 1] + alist[Nx - 1] * xlist[Nx - 2] -
             dlist[Nx - 1]) > 1.0e-6) {
        return false;
    }
}

int correctcount(bool *result, int trynumber) {
    int count = 0;
    for(int n=0; n<trynumber; n++) {
        if (result[n]) {
            count++;
        }
    }
    return count;
}

int Test(size_t diagonal_size) {
    std::chrono::time_point<std::chrono::system_clock> tstart, tend;
    std::chrono::duration<double> duration;

    PCR_Solver crs = PCR_Solver(diagonal_size);

    thrust::device_vector<float> alist(diagonal_size);
    thrust::device_vector<float> blist(diagonal_size);
    thrust::device_vector<float> clist(diagonal_size);
    thrust::device_vector<float> dlist(diagonal_size);
    thrust::device_vector<float> xlist(diagonal_size);

    float *ptr_alist = thrust::raw_pointer_cast(alist.data());
    float *ptr_blist = thrust::raw_pointer_cast(blist.data());
    float *ptr_clist = thrust::raw_pointer_cast(clist.data());
    float *ptr_dlist = thrust::raw_pointer_cast(dlist.data());
    float *ptr_xlist = thrust::raw_pointer_cast(xlist.data());

    size_t Nx = diagonal_size;
    float *a, *b, *c, *rhs, *old;
    a = new float[Nx];
    b = new float[Nx];
    c = new float[Nx];
    rhs = new float[Nx];
    old = new float[Nx];

    int trynumber = 100;
    double *gputimes = new double[trynumber];
    double *cputimes = new double[trynumber];
    bool cpuresultcheck = new bool[trynumber];
    bool gpuresultcheck = new bool[trynumber];

    double v = -1.0;
    for (int count = 0; count < trynumber; count++) {
        for (int i = 0; i < diagonal_size; i++) {
            alist[i] = -1.0 + 0.1 * float(rand()) / float(RAND_MAX);
            blist[i] = 2.0 + 0.1 * float(rand()) / float(RAND_MAX);
            clist[i] = -1.0 + 0.1 * float(rand()) / float(RAND_MAX);
            dlist[i] = 1.0 + 10.0 * float(rand()) / float(RAND_MAX);
            xlist[i] = 0.0f;
            a[i] = alist[i];
            b[i] = blist[i];
            c[i] = clist[i];
            rhs[i] = dlist[i];
            old[i] = xlist[i];
        }

        a[0] = float(0.0);
        c[Nx - 1] = float(0.0);
        alist[0] = float(0.0);
        clist[diagonal_size - 1] = float(0.0);

        tstart = std::chrono::system_clock::now();
        crs.Solve(ptr_alist, ptr_blist, ptr_clist, ptr_dlist, ptr_xlist);
        tend = std::chrono::system_clock::now();
        duration = tend - tstart;
        gputimes[count] = duration.count();
        gpuresultcheck[count] =
            ResultCheckGpu(alist, blist, clist, dlist, xlist, diagonal_size);

        tstart = std::chrono::system_clock::now();
        TDMA(a, b, c, rhs, old, Nx);
        tend = std::chrono::system_clock::now();
        duration = tend - tstart;
        cputimes[count] = duration.count();
        cpuresultcheck[count] = ResultCheckCpu(a, b, c, rhs, old, Nx);

        for (size_t it = 0; it < Nx; it++) {
            if (fabs(old[it] - xlist[it]) > v) {
                v = fabs(old[it] - xlist[it]);
            }
        }
    }

    std::cout << "Diagonal size: " << diagonal_size << std::endl;
    std::cout << "Max deviation: " << v << std::endl;
    std::cout << "gpu v.s. cpu" << std::endl;
    std::cout << "correct results: " << correctcount(gpuresultcheck, trynumber) << " <> " << correctcount(cpuresultcheck, trynumber) << std::endl;
    std::cout << "maximum: " << maxarr(gputimes, trynumber) << " <> "
              << maxarr(cputimes, trynumber) << std::endl;
    std::cout << "minimum: " << minarr(gputimes, trynumber) << " <> "
              << minarr(cputimes, trynumber) << std::endl;
    std::cout << "average: " << averarr(gputimes, trynumber) << " <> "
              << averarr(cputimes, trynumber) << std::endl;
    std::cout << std::endl;

    delete[] a;
    delete[] b;
    delete[] c;
    delete[] rhs;
    delete[] old;
    delete[] gputimes;
    delete[] cputimes;
    delete[] cpuresultcheck;
    delete[] gpuresultcheck;

    return 0;
}

int main() {
    // Generate sampel data
    srand(time(NULL));
    std::cout << "Test of the accuracy and efficiency of two algorithms that "
                 "solves the tridiagonal matrix."
              << std::endl;
    Test(10);
    Test(100);
    Test(1000);
    Test(5000);
    return 0;
}
