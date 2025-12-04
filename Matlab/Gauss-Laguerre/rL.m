function [x,w,wb]=rL(n,a,b,expoc,es)
% Gauss-Laguerre and Lobatto-Laguerre nodes and weights
% Input: 
%        n: degree
%        a: alpha
%        b: b=0 for standard Gauss and b=1 for Lobatto
% Output:
%         x: nodes
%         w: Gaussian weights
%        wb: barycentric weights (unnormalized)   
if nargin<4
  expoc=1.e3;
end    
if nargin<3
    b=0;
    if nargin<2
        a=0;
    end    
end
if (b~=0)&&(b~=1)
   disp('The third input should be 0 or 1')
   return;
end 
ap=a;
a=a+b;
[x,w]=GL(n,a,expoc);
fin=length(x);
if b==0
   si=-1;
   wb=zeros(1,fin);
   for i=1:fin
       si=-si;
       wb(i)=si*sqrt(x(i)*w(i));
   end    
else
   w=w./x;
   x=[0 x];fin=fin+1;  
   w0lo=gammaln(a)+gammaln(n+1)-gammaln(n+a+1);
   w=a*[exp(w0lo) w];
   si=1;
   wb=zeros(1,fin);
   wb(1)=sqrt(a*w(1));
   for i=2:fin
      si=-si;
      wb(i)=si*sqrt(w(i));   
   end 
end    
if(es==1)
   w=gamma(ap+1)*w;
end    
end


