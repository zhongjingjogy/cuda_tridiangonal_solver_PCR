#include <iostream>
#include <chrono>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <cstdlib>
#include <cmath>
#include "PCR_Class.h"
#include "tdma.h"

double maxarr(double *arr, int Nx) {
    double v = arr[0];
    for(int nx=0; nx<Nx; nx++) {
        if (v < arr[nx]) {
            v = arr[nx];
        }
    }
    return v;
}

double minarr(double *arr, int Nx) {
    double v = arr[0];
    for(int nx=0; nx<Nx; nx++) {
        if (v > arr[nx]) {
            v = arr[nx];
        }
    }
    return v;
}

double averarr(double *arr, int Nx) {
    double v = 0.0;
    for(int nx=0; nx<Nx; nx++) {
        v += arr[nx];
    }
    return v / double(Nx);
}

int main( ) {

    std::chrono::time_point<std::chrono::system_clock> tstart, tend;
    std::chrono::duration<double> duration;    
    size_t diagonal_size = 1000;
    
    PCR_Solver crs = PCR_Solver(diagonal_size);
    
    //Generate sampel data
    srand (time(NULL));
    
    thrust::device_vector<float> alist(diagonal_size);
    thrust::device_vector<float> blist(diagonal_size);
    thrust::device_vector<float> clist(diagonal_size);
    thrust::device_vector<float> dlist(diagonal_size);
    thrust::device_vector<float> xlist(diagonal_size);
    
    float * ptr_alist = thrust::raw_pointer_cast(alist.data());
    float * ptr_blist = thrust::raw_pointer_cast(blist.data());
    float * ptr_clist = thrust::raw_pointer_cast(clist.data());
    float * ptr_dlist = thrust::raw_pointer_cast(dlist.data());
    float * ptr_xlist = thrust::raw_pointer_cast(xlist.data());

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

    for(int count=0; count<trynumber; count++) {
        for (int i=0; i < diagonal_size; i++) {
            alist[i] = -1.0+0.1*float(rand()) / float(RAND_MAX);
            blist[i] = 2.0+0.1*float(rand()) / float(RAND_MAX);
            clist[i] = -1.0+0.1*float(rand()) / float(RAND_MAX);
            dlist[i] = 1.0 + 10.0*float(rand()) / float(RAND_MAX);
            xlist[i] = 0.0f;
            a[i] = alist[i];
            b[i] = blist[i];
            c[i] = clist[i];
            rhs[i] = dlist[i];
            old[i] = xlist[i];
        }
    
        a[0] = float(0.0);
        c[Nx-1] = float(0.0);
        alist[0] = float(0.0);
        clist[diagonal_size-1] = float(0.0);
        
        tstart = std::chrono::system_clock::now();
        crs.Solve(ptr_alist, ptr_blist, ptr_clist, ptr_dlist, ptr_xlist);
        tend = std::chrono::system_clock::now();
        duration = tend - tstart;
        gputimes[count] = duration.count();

        tstart = std::chrono::system_clock::now();
        TDMA(a, b, c, rhs, old, Nx);
        tend = std::chrono::system_clock::now();
        duration = tend - tstart;
        cputimes[count] = duration.count();

        for (size_t it=0; it<Nx; it++) {
            if(fabs(old[it] - xlist[it]) > 1.0e-6) {
                std::cout << old[it] << " " << xlist[it] << ": " << fabs(old[it] - xlist[it]) << std::endl;
            }
        }
    }

    std::cout << "maximum: " << maxarr(gputimes, trynumber) << " <> " << maxarr(cputimes, trynumber) << std::endl;
    std::cout << "minimum: " << minarr(gputimes, trynumber) << " <> " << minarr(cputimes, trynumber) << std::endl;
    std::cout << "average: " << averarr(gputimes, trynumber) << " <> " << averarr(cputimes, trynumber) << std::endl;

    delete []a;
    delete []b;
    delete []c;
    delete []rhs;
    delete []old;
    delete []gputimes;
    delete []cputimes;

    return 0;

}
