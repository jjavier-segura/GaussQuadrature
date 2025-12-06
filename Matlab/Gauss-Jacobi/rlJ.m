function [x,w,wb]=rlJ(n,a,b,e)
% Gauss-Radau-Jacobi and Gauss-Lobatto quadratures and barycentric weights
% Copyright 2025 A. Gil, J. Segura, N. M. Temme
%
% [x,w,wb]=rlJ(n,a,b,e)
%
%x: nodes, w: gauss weights, wb: barycentric weights
%n: degre, a: alpha, b: beta
%e=[c d]
%       c: set c=1 for including x=-1 as node (Gauss-Radau or Gauss-Lobatto), c=0 in other case
%       d: set d=1 for including x=1 as node (Gauss-Radau or Gauss-Lobatto), d=0 in other case
if nargin<4
   c=0;d=0;
else 
    tam=size(e);
    if (tam(1)~=1)||(tam(2)~=2)
       disp('The last argument must be a row vector with two components')
       return
    end    
    c=e(1);d=e(2);
    if ((c~=0)&&(c~=1))||((d~=0)&&(d~=1))
       disp('Wrong input in the fourth argument')
       return
    end 
end    
a=a+d;
b=b+c;
[x,w]=GJ(n,a,b);
if(c==0)&&(d==0)
   fin=length(x);
    si=-1;
    wb=zeros(1,fin);
    for i=1:fin
       si=-si;
       wb(i)=si*sqrt((1-x(i)^2)*w(i));
    end
end
if (c==1)&&(d==1)
  logw0=(a+b-1)*log(2)+gammaln(b)+gammaln(b+1)+gammaln(n+a+1)...
    +gammaln(n+1)-gammaln(n+b+1)-gammaln(n+a+b+1);
  logw1=(a+b-1)*log(2)+gammaln(a)+gammaln(a+1)+gammaln(n+b+1)...
    +gammaln(n+1)-gammaln(n+a+1)-gammaln(n+a+b+1);
  w=w./(1-x.^2);
  x=[-1 x 1];
  w=[exp(logw0) w exp(logw1)];
  wb(1)=sqrt(b*w(1));
  si=-1;
  finb=length(x);
  for i=2:finb-1
     wb(i)=si*sqrt(w(i));
     si=-si;
  end
  wb(finb)=si*sqrt(a*w(finb));
end    
if (c==1)&&(d==0)
   logw0=(a+b)*log(2)+gammaln(b)+gammaln(b+1)+gammaln(n+a+1)...
    +gammaln(n+1)-gammaln(n+b+1)-gammaln(n+a+b+1);
   w=w./(1+x);
   x=[-1 x];
   w=[exp(logw0) w];
   wb(1)=sqrt(2*b*w(1));
   si=1;
   finb=length(x);
   for i=2:finb
      si=-si; 
      wb(i)=si*sqrt((1-x(i))*w(i));
   end
end
if (c==0)&&(d==1)
    logw1=(a+b)*log(2)+gammaln(a)+gammaln(a+1)+gammaln(n+b+1)...
    +gammaln(n+1)-gammaln(n+a+1)-gammaln(n+a+b+1);
    w=w./(1-x);
    x=[x 1];
    w=[w exp(logw1)];
    si=1;
    finb=length(x);
    for i=1:finb-1 
      wb(i)=si*sqrt((1+x(i))*w(i));
      si=-si;
    end
    wb(finb)=si*sqrt(2*a*w(finb));
end
