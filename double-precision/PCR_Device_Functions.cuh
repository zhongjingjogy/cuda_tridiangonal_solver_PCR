#ifndef PCR_Device_Functions_cuh
#define PCR_Device_Functions_cuh

__global__ void list_print(int nmax, double *in);

__global__ void Solve_Kernel(double * alist, double * blist, double * clist, double * dlist, double * xlist, int iter_max, int DMax);


#endif
