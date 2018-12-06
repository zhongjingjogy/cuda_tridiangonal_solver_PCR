# Tridiagonal matrix equation
Tridiagonal matrix equation is one very important equation for numerical simulation, where the implicit finite difference is applied to solve PDE.

# cuda_tridiangonal_solver_PCR
Implement Parallel Cyclic Reduction in CUDA-C for a Tridiagonal matrix equation

# Measurement of performance
A traditional TDMA algorithm is implemented in order to compare with the current PCR tridiagonal solver.

The test environment is `2 vCPU 8 GB, GPU：Nvidia Tesla P4`. The two algorithms are applied to solve a linear equation `Ax=d`, where `A` is a diagonal dominant matrix. The matrix `A` and the right hand side are both generated randomly. The two algorithms will run repeatedly, and the statistics about the elapsed time will be used as the measurement of the efficiency. Meanwhile, comparison between the calculated results are made directly as a measurement of the accuracy.

Here comes the results,

```
Test of the accuracy and efficiency of two algorithms that solves the tridiagonal matrix.
Diagonal size: 10
Max deviation: 1.90735e-05
gpu v.s. cpu
maximum: 1.0491e-05 <> 4.193e-06
minimum: 4.593e-06 <> 3.59e-07
average: 5.84073e-06 <> 5.3118e-07

Diagonal size: 100
Max deviation: 5.34058e-05
gpu v.s. cpu
maximum: 9.163e-06 <> 3.61e-06
minimum: 5.37e-06 <> 2.65e-06
average: 7.19839e-06 <> 2.88266e-06

Diagonal size: 1000
Max deviation: 4.95911e-05
gpu v.s. cpu
maximum: 1.0157e-05 <> 4.1617e-05
minimum: 7.829e-06 <> 2.1209e-05
average: 9.12968e-06 <> 2.20262e-05

Diagonal size: 5000
Max deviation: 87.031
gpu v.s. cpu
maximum: 5.765e-06 <> 0.000113033
minimum: 2.044e-06 <> 0.000101797
average: 2.33923e-06 <> 0.000103271
```

It might be found that as the size of the diagonal matrix increases, the efficiency of the one with GPU implementation goes above that with CPU version. But it seems that the deviation between the two algorithms are not negligible. The detail implementation should be checked for both the algorithms.