function [xc,w]=bH(n,expoc)
% Hermite nodes and barycentric weights. Iterative method.
% Copyright 2025 A. Gil, J. Segura, N. M. Temme
%
% [x,v]=bH(n,expoc)
%
% INPUT:
%   n: degree 
%   expoc (optional): cut off (only weights such that v/vmax<10^(-expoc) are computed)
%
% OUTPUT:
%   x: nodes
%   v: barycentric weights
%------------------------------------------------------------
% This algorithm computes barycentric weights directly in
% terms of the derivative of the OP instead of relating them
% to Gauss weights. The range of computation is extended.
% The normalization is different from that in GHi.m (the 
% normalization is arbitrary in the barycentric formula)
%------------------------------------------------------------
if nargin<2 
    expoc=1.e300;
end
pi=3.1415926535897932385;
epsil=1.0e-9;
Nn=floor(n/2);
ind=n-2*Nn;
ta=floor(n/2.d0)+ind;
xc=zeros(1,ta);deri=zeros(1,ta);w=zeros(1,ta);
xp=0;
delta=(1+ind)/2*pi/sqrt(2*n+1);
xn=xp+delta;
if ind>0
  u0=0;
  up1=1;
  i=1;
  xc(i)=0;
  deri(i)=1;
else
  u0=1;
  up1=0;
  i=0;
end
while i<Nn+ind && xn<sqrt(4.6*expoc)
  i=i+1;
  err=1+epsil;
  while err>epsil
    t=xn-xp;
    delta=t;
    um1=0;
    um2=0;
    suma=u0+up1*t;
    sumad=up1;
    k=-1;
    errod=1;
    while ((errod>1.0e-25)&&(k<50))||(k<5)
      k=k+1;
      up2=(xp*xp-2*n-1)*u0+2*k*xp*um1+k*(k-1)*um2;
      um2=um1;
      um1=u0;u0=up1;up1=up2;
      ntd=up2*t;
      sumad=sumad+ntd;
      if (i>1+ind) 
        errod=abs(ntd/sumad);
      end    
      t=t*delta/(k+2);
      nt=up2*t;
      suma=suma+nt;
    end
    u0=suma;
    up1=sumad;
    hs=u0/up1;
    ws=sqrt(2*n+1-xn*xn);
    delta=-atan(ws*hs)/ws;
    xp=xn;
    xn=xp+delta;       
    err=abs(delta/xn);
  end
  xc(i)=xn;
  t=xn-xp;
  delta=t;
  um1=0;
  um2=0;
  suma=u0+up1*t;
  sumad=up1;
  k=-1;
  errod=1;
  while ((errod>1.0e-25)&&(k<50))||(k<5)
    k=k+1;
    up2=(xp*xp-2*n-1)*u0+2*k*xp*um1+k*(k-1)*um2;
    um2=um1;
    um1=u0;u0=up1;up1=up2;
    ntd=up2*t;
    sumad=sumad+ntd;
    if i>1+ind 
      errod=abs(ntd/sumad);
    end    
    t=t*delta/(k+2);
    nt=up2*t;
    suma=suma+nt;
  end
  u0=suma;
  up1=sumad;
  deri(i)=sumad;
  delta=pi/sqrt(2*n+1-xn*xn);
  xp=xn;
  xn=xp+delta;
end
longi=i;
for j=1:longi
  w(j)=exp(-0.5*xc(j)^2)/deri(j);
end
if n/2==floor(n/2)
   xc=[-xc(longi:-1:1) xc(1:longi)];  
   w=[-w(longi:-1:1) w(1:longi)];
else    
   xc=[-xc(longi:-1:2) xc(1:longi)];  
   w=[w(longi:-1:2) w(1:longi)];
end   
