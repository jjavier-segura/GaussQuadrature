function [xc,wns,w,uv]=GHi(n,expoc)
%-----------------------------------------------------------------------
%Iterative algorithm.
%Double precision Hermite nodes and weights
% n-->degree 
% expoc: cut off (only weights such that w/wmax<10^(-expoc) are computed
% xc-->nodes;  wns-->unscaled weights; 
% w-->scaled weights; uv-->barycentric weights
%-----------------------------------------------------------------------
if nargin<2 
    expoc=1.e300;
end
expoc=max(expoc,20);
pi=3.1415926535897932385;
epsil=1.0e-9;
Nn=floor(n/2);
ind=n-2*Nn;
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
xcut=sqrt(2.3*expoc);
while i<Nn+ind && xn<xcut
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
istop=min(i,Nn+ind);
sumap=0;
i=0;
xco=0;
while (xco<20)&&(i<istop)
  i=i+1;
  xp2=xc(i)^2;
  xco=xc(i);
  sumap=sumap+2*xp2*exp(-xp2)/deri(i)^2;
end
i=1:istop;
w=0.5*sqrt(pi)/sumap./deri(i).^2;
wns=w.*exp(-xc(i).^2);
xc=[-xc(end:-1:1+ind) xc];
wns=[wns(end:-1:1+ind) wns];
w=[w(end:-1:1+ind) w];
uv=-cumprod(-ones(1,length(xc))).*sqrt(wns);
