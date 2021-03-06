#include "PCR_Class.h"
#include "PCR_Device_Functions.cuh"


int count_iter(int bin) {

	int m = 0;
	int count = 0;

	while (bin != 0) {
		m+=bin%2;
		bin = bin >> 1;
		count++;
	}
	count--;

	if (m>1) m = 0;

    return (count - m);
}

PCR_Solver::PCR_Solver(int coming_ds) {

    diagonal_size = coming_ds;

    iter_max = count_iter(coming_ds);

}


void PCR_Solver::Solve(double * alist, double * blist, double * clist, double * dlist, double * xlist) {

	Solve_Kernel<<<1,diagonal_size>>>(alist, blist, clist, dlist, xlist, iter_max, diagonal_size);

}
