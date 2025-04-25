quadg:=proc(f,n,lambda)
local N,ind,i,In:
N:=floor(n/2):ind:=n-2*N:
ggnw(n,lambda,x,w):
In:=0:
for i from 1-ind to N do:
 In:=In+evalf(w[i]*subs(x=x[i],f)):
end:
for i from 1 to N do:
 In:=In+evalf(w[i]*subs(x=-x[i],f)):
end:

end proc:


ggnw:=proc(n,lambda,xc,w)
local epsil,N,ind,epso,L,xp,delta,xn,y0,yp1,err,ym1,
ym2,t,suma,sumad,j,errod,c1,c2,c3,c4,c5,yp2,ntd,deri,ws,sumws,i,
nt,hs,de,sumw,kk,factor;
description "Gauss-Gegenbauer quadrature";
epsil:=evalf(10^(-Digits+4)):
N:=floor(n/2):ind:=n-2*N:epso:=evalf(10^(-Digits)):
L:=2*(n+lambda)+1:
xp:=0:
delta:=evalf(tanh((1+ind)*Pi/sqrt(L^2-4*lambda^2-1))):
xn:=evalf(xp+delta):
if ind>0 then
   y0:=0;yp1:=1;
else
   y0:=1;yp1:=0;
end if:
i:=0:   
while i<N do
    i:=i+1:
    err:=1+epsil:
    while err>epsil do
      t:=delta;
      ym1:=0;ym2:=0; 
      suma:=evalf(y0+yp1*t):
      sumad:=yp1;
      j:=-1;
      errod:=1; 
      while errod>epso and j<4000 do
          j:=j+1;
          c1:=evalf(4*(j+2)*(j+1)*(1-xp^2)^2):
          c2:=evalf(-16*(j+1)*j*(1-xp^2)*xp):
          c3:=evalf(0.5*j*(j-1)*(48*xp^2-16)+(L*L-1)*(1-xp^2)+4*(1-lambda^2)):
          c4:=evalf(16*(j-1)*(j-2)*xp-2*(L*L-1)*xp):
          c5:=evalf(4*(j-2)*(j-3)+(1-L*L)):
          yp2:=evalf(-1/c1*(c2*yp1+c3*y0+c4*ym1+c5*ym2)):
          ym2:=ym1;ym1:=y0;y0:=yp1;yp1:=yp2:  
          ntd:=(j+2)*yp2*t;
          sumad:=evalf(sumad+ntd);
          t:=t*delta;
          nt:=yp2*t;
          suma:=evalf(suma+nt);
          if i>1 and j>20 then
             errod:=max(abs(evalf(ntd/sumad)),abs(evalf(nt/suma)));
          end;
      end do;
      y0:=suma;yp1:=sumad:
      hs:=y0/(yp1*(1-xn^2)+xn*y0):
      ws:=evalf(0.5*sqrt((1-L^2)*xn^2+L^2-1-4*lambda^2));
      de:=evalf(tanh(arctan(ws*hs)/ws));
      xp:=xn;
      xn:=(xp-de)/(1-xp*de);     
      delta:=xn-xp:  
      err:=abs(delta/xn);
    end do;
    xc[i]:=xn:
    t:=delta;
    ym1:=0;ym2:=0; 
    suma:=evalf(y0+yp1*t):
    sumad:=yp1;
    j:=-1:
    errod:=1; 
    while errod>epso and j<4000 do
           j:=j+1;
           c1:=evalf(4*(j+2)*(j+1)*(1-xp^2)^2):
           c2:=evalf(-16*(j+1)*j*(1-xp^2)*xp):
           c3:=evalf(0.5*j*(j-1)*(48*xp^2-16)+(L*L-1)*(1-xp^2)+4*(1-lambda^2)):
           c4:=evalf(16*(j-1)*(j-2)*xp-2*(L*L-1)*xp):
           c5:=evalf(4*(j-2)*(j-3)+(1-L*L)):
          yp2:=evalf(-1/c1*(c2*yp1+c3*y0+c4*ym1+c5*ym2));
          ym2:=ym1;ym1:=y0;y0:=yp1;yp1:=yp2:  
          ntd:=(j+2)*yp2*t;
          sumad:=evalf(sumad+ntd);
          t:=t*delta;
          nt:=yp2*t;
          suma:=evalf(suma+nt);
          if i>1 and j>20 then
             errod:=max(abs(evalf(ntd/sumad)),abs(evalf(nt/suma)));
          end;
    end do:
    y0:=suma;yp1:=sumad;
    deri[i]:=sumad;
    de:=evalf(tanh(Pi/ws));
    xp:=xn;
    xn:=(xp+de)/(1+xp*de);     
    delta:=xn-xp:
end do:
sumw:=0:
for i from 1 to N do
     ws[i]:=evalf((1-xc[i]^2)^(lambda)/deri[i]^2);
     sumw:=evalf(sumw+ws[i]);
end:
if ind=1 then
   xc[0]:=0;
   deri[0]:=1;
   ws[0]:=1;
   sumw:=evalf(sumw+ws[0]/2);
end:
factor:=evalf(0.5*sqrt(Pi))*evalf(GAMMA(lambda + 1)/GAMMA(lambda + 3/2));
for i from 1-ind to N do
     w[i]:=factor*ws[i]/sumw:
end:
end proc;
