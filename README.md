# Tridiagonal matrix equation
Tridiagonal matrix equation is one very useful equation for numerical simulation, where the implicit finite difference is applied to solve PDE. This repository is to go through different implementations of the TDMA algorithms and try to figure out the efficiency and accuracy of the algorithms.

## Accuracy test on the CPU implementation
The single-precision and double-precision float number are used to store the coefficients of the tridiagonal test. Then the tridiagonal matrix algorithm will be used to solve the equation. Two implementations of the TDMA are presented here, i.e., Wiki and mine. 

Test environment: MinGW-W64, Windows 10
```
g++ (x86_64-posix-sjlj-rev0, Built by MinGW-W64 project) 6.4.0
Copyright (C) 2017 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
```

Here come the test result,
```
Test of the accuracy and efficiency of tdma algorithm that solves the tridiagonal matrix.
Test 1....
Max absolute error: 9.53674e-007
Elapsed: 0 second
0.148776 0.75612 -1.00188 2.25141 

Test 2: accuracy test on single-precision float....
size of a, b, c, rhs, x: 100 100 100 100 100
Wiki implementation: ...
Comparison with the standard result: 7 false 21 false 53 false 58 false 59 false 60 false 86 false 98 false 
Max absolute error: 3.8147e-006
Jing's implementation: ...
Comparison with the standard result: 7 false 21 false 53 false 58 false 59 false 60 false 86 false 98 false 
Max absolute error: 3.8147e-006

Test 3: accuracy test on double-precision float....
size of a, b, c, rhs, x: 100 100 100 100 100
Wiki implementation: ...
Comparison with the standard result: 
Max absolute error: 8.88178e-015
Jing's implementation: ...
Comparison with the standard result: 
Max absolute error: 8.88178e-015

Test 4: accuracy test on double-precision float....
Diagonal size: 10
correct results: 100 / 100, tolerance: 1e-010

Test 5: accuracy test on double-precision float....
Diagonal size: 100
correct results: 100 / 100, tolerance: 1e-010

Test 6: accuracy test on double-precision float....
Diagonal size: 500
correct results: 100 / 100, tolerance: 1e-010

Test 7: accuracy test on double-precision float....
Diagonal size: 1000
correct results: 100 / 100, tolerance: 1e-010

Test 8: accuracy test on double-precision float....
Diagonal size: 10000
correct results: 100 / 100, tolerance: 1e-010

Test 9: accuracy test on single-precision float....
Diagonal size: 10
correct results: 25 / 100, tolerance: 1e-005

Test 10: accuracy test on single-precision float....
Diagonal size: 100
correct results: 0 / 100, tolerance: 1e-005

Test 11: accuracy test on single-precision float....
Diagonal size: 500
correct results: 0 / 100, tolerance: 1e-005

Test 12: accuracy test on single-precision float....
Diagonal size: 1000
correct results: 0 / 100, tolerance: 1e-005

Test 13: accuracy test on single-precision float....
Diagonal size: 10000
correct results: 0 / 100, tolerance: 1e-005

```

## Efficiency and accuracy of the Tridiagonal matrix equation solved by Parallel Cyclic Reduction

A traditional TDMA algorithm is implemented in order to compare with the current PCR tridiagonal solver.

The test environment is `2 vCPU 8 GB, GPU：Nvidia Tesla P4`. The two algorithms are applied to solve a linear equation `Ax=d`, where `A` is a diagonal dominant matrix. The matrix `A` and the right hand side are both generated randomly. The two algorithms will run repeatedly, and the statistics about the elapsed time will be used as the measurement of the efficiency. Meanwhile, comparison between the calculated results are made directly as a measurement of the accuracy.

Here comes the results,

```
Test of the accuracy and efficiency of two algorithms that solves the tridiagonal matrix.
Diagonal size: 10
Gpu v.s. Cpu
Correct results: 100 <> 100
Maximum time(s): 1.081e-05 <> 4.407e-06
Minimum time(s): 3.955e-06 <> 4.36e-07
Average time(s): 5.76189e-06 <> 6.4576e-07

Diagonal size: 100
Gpu v.s. Cpu
Correct results: 100 <> 100
Maximum time(s): 1.8211e-05 <> 3.427e-06
Minimum time(s): 5.493e-06 <> 2.857e-06
Average time(s): 7.28943e-06 <> 3.1191e-06

Diagonal size: 1000
Gpu v.s. Cpu
Correct results: 58 <> 100
Maximum time(s): 1.4395e-05 <> 3.3783e-05
Minimum time(s): 7.772e-06 <> 2.3135e-05
Average time(s): 8.99534e-06 <> 2.35596e-05

Diagonal size: 5000
Gpu v.s. Cpu
Correct results: 0 <> 100
Maximum time(s): 6.806e-06 <> 0.000154798
Minimum time(s): 2.04e-06 <> 0.000112074
Average time(s): 2.25151e-06 <> 0.000113564
```

It might be found that as the size of the diagonal matrix increases, the efficiency of the one with GPU implementation goes above that with CPU version. But it seems that the deviation between the two algorithms are not negligible. The detail implementation should be checked for both the algorithms.