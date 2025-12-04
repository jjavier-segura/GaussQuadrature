# Gauss quadrature
This project provides fast algorithms for computing Gauss-Jacobi (Gauss-Legendre as a particular case), Gauss-Laguerre and Gauss-Hermite quadratures together with Gauss-Radau and Gauss-Lobatto variants. In addition to Matlab algorithms, arbitrary precision Maple algorithms are provided for the symmetric cases (Gauss-Gegenbauer and Hermite). 

The algorithms are described in the paper "Fast and accurate computation of classical Gaussian quadratures", by A. Gil, J. Segura & N. M. Temme (submitted to SIAM J. Sci. Comput.)

![image](https://github.com/user-attachments/assets/bd4bd13d-e3c0-48fc-9f11-49f7df991eda)

# Overview of the software

# Gauss-Jacobi quadratures

```matlab
function [x,w,ie]=GJ(n,a,b,me)
```
Computation of Gauss-Jacobi quadrature (uses the function `legen`). 

```matlab
function [x,w]=legen(n)
```
Fast computation of Gauss-Legendre quadrature.

```matlab
function [x,w,wb]=rlJ(n,a,b,e)
```
Gauss-Radau-Jacobi and Gauss-Lobatto-Jacobi quadratures and barycentric weights (uses the function `GJ`). 


