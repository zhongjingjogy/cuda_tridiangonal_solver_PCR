#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "pcr/PCR_Class.h"
#include "tdma/tdma.h"
#include "utils.h"

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
            a[i] = alist[i];
            b[i] = blist[i];
            c[i] = clist[i];
            rhs[i] = dlist[i];
            old[i] = xlist[i];
        }

        a[0] = double(0.0);
        c[Nx - 1] = double(0.0);
        alist[0] = float(0.0);
        clist[diagonal_size - 1] = float(0.0);

        tstart = std::chrono::system_clock::now();
        crs.Solve(ptr_alist, ptr_blist, ptr_clist, ptr_dlist, ptr_xlist);
        tend = std::chrono::system_clock::now();
        duration = tend - tstart;
        gputimes[count] = duration.count();
        gpuresultcheck[count] =
            ResultCheckGpu(alist, blist, clist, dlist, xlist, diagonal_size) < 1.0e-5;

        tstart = std::chrono::system_clock::now();
        TDMA(a, b, c, rhs, old, Nx);
        tend = std::chrono::system_clock::now();
        duration = tend - tstart;
        cputimes[count] = duration.count();
        cpuresultcheck[count] = ResultCheckCpu(a, b, c, rhs, old, Nx) < 1.0e-15;
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
