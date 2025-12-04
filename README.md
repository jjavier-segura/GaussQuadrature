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

# Gauss-Laguerre quadratures

```matlab
[x,w,we,indi,ie]=GL(n,a,expoc,es,me)
```
Computation of Gauss-Laguerre quadrature (uses the function `bessJY`).

```matlab
function [x,w,wb]=rL(n,a,b,expoc)
```
Computation of Gauss-Radau-Laguerre quadrature and barycentric weights (uses the function `GL`).

# Gauss-Hermite quadrature

```matlab
function [xc,wns,w,uv]=GHa(n,expoc)
```
Asymptotic computation of Gauss-Hermite quadrature.

```matlab
function [xc,wns,w,uv]=GHi(n,expoc)
```
Iterative computation of Gauss-Hermite quadrature.

```matlab
function [xc,wb]=bH(n,expoc)
```

Computation of barycentric weights. 
This algorithm computes barycentric weights directly in terms of the derivative of the orthogonal polynomial instead of relating them to Gauss weights. The range of computation is extended.
 





