#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "PCR_Class.h"
#include "tdma.h"
#include "utils.h"

float GpuResultCheck(float *a, float *b, float *c, float *rhs, float *old, int Nx) {
    float v = -1.0;
    float v1 = fabs(b[0] * old[0] + c[0] * old[1] - rhs[0]);
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
    double *a, *b, *c, *rhs, *old;
    a = new double[Nx];
    b = new double[Nx];
    c = new double[Nx];
    rhs = new double[Nx];
    old = new double[Nx];

    float *ma, *mb, *mc, *mrhs, *mx;
    ma = new float[Nx];
    mb = new float[Nx];
    mc = new float[Nx];
    mrhs = new float[Nx];
    mx = new float[Nx];

    int trynumber = 100;
    double *gputimes = new double[trynumber];
    double *cputimes = new double[trynumber];
    bool *cpuresultcheck = new bool[trynumber];
    bool *gpuresultcheck = new bool[trynumber];

    for (int count = 0; count < trynumber; count++) {
        for (int i = 0; i < diagonal_size; i++) {
            alist[i] = -1.0 + 0.1 * float(rand()) / float(RAND_MAX);
            blist[i] = 2.0 + 0.1 * float(rand()) / float(RAND_MAX);
            clist[i] = -1.0 + 0.1 * float(rand()) / float(RAND_MAX);
            dlist[i] = 1.0 + 10.0 * float(rand()) / float(RAND_MAX);
            xlist[i] = 0.0f;

            ma[i] = alist[i];
            mb[i] = blist[i];
            mc[i] = clist[i];
            mrhs[i] = dlist[i];
            mx[i] = xlist[i];

            a[i] = double(ma[i]);
            b[i] = double(mb[i]);
            c[i] = double(mc[i]);
            rhs[i] = double(mrhs[i]);
            old[i] = double(mx[i]);
        }

        a[0] = double(0.0);
        c[Nx - 1] = double(0.0);
        alist[0] = float(0.0);
        clist[diagonal_size - 1] = float(0.0);
        ma[0] = float(0.0);
        mc[diagonal_size-1] = float(0.0);

        tstart = std::chrono::system_clock::now();
        crs.Solve(ptr_alist, ptr_blist, ptr_clist, ptr_dlist, ptr_xlist);
        tend = std::chrono::system_clock::now();
        duration = tend - tstart;
        gputimes[count] = duration.count();
        for (int i = 0; i < diagonal_size; i++) {
            mx[i] = xlist[i];
        }
        gpuresultcheck[count] =
            GpuResultCheck(ma, mb, mc, mrhs, mx, diagonal_size) < 1.0e-4;

        tstart = std::chrono::system_clock::now();
        TDMA<double>(a, b, c, rhs, old, Nx);
        tend = std::chrono::system_clock::now();
        duration = tend - tstart;
        cputimes[count] = duration.count();
        cpuresultcheck[count] = CpuResultCheck(a, b, c, rhs, old, Nx) < 1.0e-10;
    }

    std::cout << "Diagonal size: " << diagonal_size << std::endl;
    std::cout << "Gpu v.s. Cpu" << std::endl;
    std::cout << "Correct results: " << correctcount(gpuresultcheck, trynumber) << " <> " << correctcount(cpuresultcheck, trynumber) << std::endl;
    std::cout << "Maximum time(s): " << maxarr(gputimes, trynumber) << " <> "
              << maxarr(cputimes, trynumber) << std::endl;
    std::cout << "Minimum time(s): " << minarr(gputimes, trynumber) << " <> "
              << minarr(cputimes, trynumber) << std::endl;
    std::cout << "Average time(s): " << averarr(gputimes, trynumber) << " <> "
              << averarr(cputimes, trynumber) << std::endl;
    std::cout << std::endl;

    delete[] a;
    delete[] b;
    delete[] c;
    delete[] rhs;
    delete[] old;
    delete[] ma;
    delete[] mb;
    delete[] mc;
    delete[] mrhs;
    delete[] mx;
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
