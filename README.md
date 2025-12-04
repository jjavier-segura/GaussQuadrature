# GAUSS QUADRATURE

This project provides fast algorithms for computing Gauss-Jacobi (Gauss-Legendre as a particular case), Gauss-Laguerre and Gauss-Hermite quadratures together with Gauss-Radau and Gauss-Lobatto variants. In addition to Matlab algorithms, arbitrary precision Maple algorithms are provided for the symmetric cases (Gauss-Gegenbauer and Hermite). 

The algorithms are described in the paper "Fast and accurate computation of classical Gaussian quadratures", by A. Gil, J. Segura & N. M. Temme (submitted to SIAM J. Sci. Comput.)
------------------------------------------------------------------------------------------
![image](https://github.com/user-attachments/assets/bd4bd13d-e3c0-48fc-9f11-49f7df991eda)
------------------------------------------------------------------------------------------

# Overview of the software

# MATLAB algorithms

## 1. Gauss--Jacobi quadratures

### 1.1. `GJ`

``` matlab
function [x,w,ie] = GJ(n,a,b,me)
```

Computation of Gauss--Jacobi quadrature (uses `legen`).

**Inputs:** - `n` (degree) - `a`, `b` (parameters of the quadrature) 
-`me` (mode of computation, default `0`, combined method)

**Outputs:**  
- `x` (nodes) - `w` (weights) 
- `ie` (error flag): `0` computation succesful
                     `1` computation failed
                     `2` some weights are too small and are set to zero
-----------------------------------------------------------------------------------

### 1.2. `legen`

``` matlab
function [x,w] = legen(n)
```

Fast computation of Gauss--Legendre quadrature.

**Input:**
- `n` (degree)

**Outputs:**
- `x` (nodes)
- `w` (weights)

------------------------------------------------------------------------

### 1.3. `rlJ`

``` matlab
function [x,w,wb] = rlJ(n,a,b,e)
```

Gauss--Radau--Jacobi and Gauss--Lobatto--Jacobi quadratures and
barycentric weights (uses `GJ`).

**Inputs:** - `n` (degree) - `a`, `b` (parameters) - `e = [c d]`: -
`c = 1`: include `x = -1` - `c = 0`: otherwise - `d = 1`: include
`x = 1` - `d = 0`: otherwise

**Outputs:** - `x` (nodes) - `w` (weights) - `wb` (barycentric weights)

------------------------------------------------------------------------

## 2. Gauss--Laguerre quadratures

### 2.1. `GL`

``` matlab
[x,w,we,indi,ie] = GL(n,a,expoc,es,me)
```

Computation of Gauss--Laguerre quadrature (uses `bessJY`).

**Inputs:** - `n` (degree) - `a` (parameter) - `expoc` (exponential
cutoff) - `es` (extra scaling) - `me` (optional, default `0`)

**Outputs:** - `x` (nodes) - `w` (weights) - `we` (scaled weights) -
`indi` (pivotal index) - `ie` (error flag)

------------------------------------------------------------------------

### 2.2. `rL`

``` matlab
function [x,w,wb] = rL(n,a,b,expoc)
```

Gauss--Radau--Laguerre quadrature and barycentric weights (uses `GL`).

**Inputs:** - `n` (degree) - `a` (parameter) - `b` (Gauss or Lobatto
option) - `expoc` (cutoff)

**Outputs:** - `x` (nodes) - `w` (weights) - `wb` (barycentric weights,
unscaled)

------------------------------------------------------------------------

## 3. Gauss--Hermite quadrature

**Common Inputs:**
- `n` (degree)
- `expoc` (cutoff)

### 3.1. `GHa` and `GHi`

``` matlab
function [xc,wns,w,uv] = GHa(n,expoc)
function [xc,wns,w,uv] = GHi(n,expoc)
```

Asymptotic (`GHa`) and iterative (`GHi`) computation.

**Outputs:** - `xc` (nodes) - `wns` (unscaled weights) - `w` (scaled
weights) - `uv` (barycentric weights)

------------------------------------------------------------------------

### 3.2. `bH`

``` matlab
function [xc,wb] = bH(n,expoc)
```

Direct computation of barycentric weights.

**Outputs:** - `xc` (nodes) - `wb` (barycentric weights)

------------------------------------------------------------------------

# MAPLE algorithms

1.  **Gauss--Gegenbauer quadrature:**
    -   `gegenbauer.mpl` (Maple algorithm)
    -   `gegenbauer.mws` (usage worksheet)
2.  **Gauss--Hermite quadrature:**
    -   `hermite.mpl` (Maple algorithm)
    -   `hermite.mws` (usage worksheet)


 



