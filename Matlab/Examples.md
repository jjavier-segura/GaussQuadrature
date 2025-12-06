**EXAMPLE 1**

``` matlab
%Example script #1
% Use of the function for the Gauss-Legendre quadrature 
format long e
n = 10000; % Set n=10000
[x, w] = legen(n);  % Generate the Legendre nodes and weights
disp('First positive node')
x(5001)  % Display the first positive node
disp('and corresponding weight')
w(5001)  % Display the corresponding weight for the node
disp('Check: sum of the weights')
sum(w)  % Check the sum of the weights

% --------Expected results---------------------
%First positive node
% 1.570717782483478e-04
%
%and corresponding weight
% 3.141435539132268e-04
%
%Check: sum of the weights
% 1.999999999999999e+00
```

**EXAMPLE 2**

``` matlab

% Example script #2
% Use of the function for the Gauss-Jacobi quadrature
format long e
% Degree of the quadrature
n = 1000;
% Parameters of the quadrature
a = 2; 
b = 50;
[x, w] = GJ(n, a, b);  % Generate the Jacobi nodes and weights
                        % Default computation

disp('Last computed node')
x(end)  % Display the last computed node
disp('Minimum weight computed')
min(w)  % Display the minimum weight computed
disp('Maximum weight computed')
max(w)  % Display the maximum weight computed

%--------Expected results---------------------
% Last computed node
% 9.999874773822709e-01
%
% Minimum weight computed
% 4.280602158063998e-144
%
% Maximum weight computed
% 1.161328831340990e+09
```
**EXAMPLE 3**

``` matlab

% Example script #3
% Use of the function for the Gauss-Jacobi quadrature
format long e
% Degree of the quadrature
n = 10000;
% Parameters of the quadrature
a = 90; 
b = -0.5;
[x, w, ie] = GJ(n, a, b);  % Generate the Jacobi nodes and weights

disp('Error flag')
ie %Error flag
disp('Last computed node')
x(end)  % Display the last computed node
disp('Minimum weight computed')
min(w)  % Display the minimum weight computed
disp('Maximum weight computed')
max(w)  % Display the maximum weight computed

%----------------Expected results---------------------------
% Warning: some weights are too small and have been set to zero
%
% Error flag
% ie = 2
%
% Last computed node
% 9.999518791384308e-01
%
% Minimum weight computed
% 0
%
% Maximum weight computed
% 5.475299093998360e+23
```
**EXAMPLE 4**

``` matlab

% Example script #4
% Use of the function for the Gauss-Laguerre quadrature
format long e
% Degree of the quadrature
n = 100;
% Parameter of the quadrature
a = 200; 
[x, w] = GL(n, a);  % Generate the Laguerre nodes and weights
                    % Default computation

disp('First computed node')
x(1)  % Display the first computed node
disp('Last computed node')
x(end)  % Display the last computed node
disp('Sum of weights')
sum(w)  % Display the sum of the weights
%---------------------------------------------------------
% Note that the sum of the weights is normalized to one.
% For the standard normalization, multiply the weights by
% Gamma(a+1) or choose es=1 (see Example 5).

% ----------Expected results----------------------------
% First computed node
%  5.966112947657201e+01
%
% Last computed node
%  7.143181908293859e+02
%
% Sum of weights
%  1
```
**EXAMPLE 5**

``` matlab

% Example script #5
% Use of the function for the Gauss-Laguerre quadrature
format long e
% Degree of the quadrature
n = 1000;
% Parameter of the quadrature
a = 100; 
% Exponential cutoff for the weights
expoc = 100;
% Parameter for scaling:
%   es = 0, weights normalized to 1
%   es = 1, standard normalization 
%         (normalized to Gamma(a+1))
es = 1;
[x, w, we] = GL(n, a, expoc, es);  % Laguerre nodes and weights

disp('Number of calculated nodes (or weights)')
length(x)  % Display the number of nodes (or weights) calculated
% Since expoc = 100, only the weights (and corresponding nodes)
% such that, approximately, w/max(w) < 10^(-expoc) are computed

disp('Sum of weights')
sum(w)  % Display the sum of the weights
% Note that now the sum of the weights is normalized to Gamma(a+1)

disp('Minimum of weights / Maximum of weights')
min(w) / max(w)  % Display the ratio of the minimum to maximum weight
% This value should be close to 10^(-expoc)

% ----------Expected results----------------------------
% Number of calculated nodes (or weights)
% 396
%
% Sum of weights
%  9.332621544394759e+157
%
% Minimum of weights / Maximum of weights
%  4.099853017526646e-101
```
**EXAMPLE 6**

``` matlab

% Example script #6
% Use of the function for the Gauss-Laguerre quadrature
%---------------------------------------------
% In this example, we show the role of the
% different scaling options
%---------------------------------------------
format long e
% Degree of the quadrature
n = 1000;
% Parameter of the quadrature
a = 1000; 
% Exponential cutoff for the weights
expoc = 300;
% Parameter for scaling:
%   es = 0, weights normalized to 1
%   es = 1, standard normalization 
%         (normalized to Gamma(a+1))
es = 1;

disp('Case 1')
[x, w] = GL(n, a, expoc, es);  
% An error occurs because Gamma(a+1) overflows.

% Choosing es = 0 solves the problem.
es = 0;
disp('Case 2')
[x, w] = GL(n, a, expoc, es);  
% Not all nodes and weights are computed
length(x)  % Display the number of nodes (or weights) calculated

disp('Case 3')
% If expoc is not used as an explicit input, no cutoff is applied
[x, w, we, indi] = GL(n, a);  
% A warning appears because some weights are too small
disp('Minimum of the weights')
min(w)
disp('Minimum of the scaled weights')
min(we)
disp('Maximum of the scaled weights')
max(we)
disp('Index for the pivotal node used in scaling')
indi

% ----------Expected results----------------------------
% Case 1
% Error. Degree or parameters out of range
% es = 0
%
% Case 2
% 687
%
% Case 3
% Warning: some weights are too small and have been set to zero
% Minimum of the weights
%  0
% Minimum of the scaled weights
%  3.962026614996834e-02
% Maximum of the scaled weights
%  3.069980359731089e-01
% Index for the pivotal node used in scaling
%  319
```
**EXAMPLE 7**

``` matlab

% Example 7
% Use of the function for the Gauss-Hermite quadrature 
% The function implementing the iterative method is used (function GHi)
format long e
n = 10000; % Set n = 10000
expoc = 100; % Cut-off
[x, w, ws] = GHi(n, expoc);  % Calculate the Hermite nodes and weights
                              % (unscaled and scaled)
disp('Number of nodes (or corresponding weights) calculated')
length(x)  % Display the number of nodes (or corresponding weights) calculated
disp('Minimum of the weights / Maximum of the weights')
min(w) / max(w)  % Display the ratio of the minimum to maximum weight
disp('Weight number 100, w(100)')
w(100)  % Display weight number 100
disp('Scaled weight number 100, ws(100)')
ws(100)  % Display scaled weight number 100
disp('The relation between w(100) and ws(100) is given by')
disp('ws(100)*exp(-x(100)^2)')
disp('Check:')
ws(100) * exp(-x(100)^2)  % Verify the relationship between w(100) and ws(100)

% --------Expected results---------------------
% Number of nodes (or corresponding weights) calculated
% 1362
%
% Minimum of the weights / Maximum of the weights
% 2.404346469225838e-100
%
% Weight number 100, w(100)
% 4.787503794974321e-75
%
% Scaled Weight number 100, ws(100)
% 2.230736516171637e-02
%
% The relation between w(100) and ws(100) is given by
% ws(100)*exp(-x(100)^2)
% Check:
% 4.787503794974321e-75
```
