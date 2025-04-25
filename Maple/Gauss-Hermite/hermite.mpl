quadh:=proc(f,n)
local N,ind,i,In:
N:=floor(n/2):ind:=n-2*N:
ghnw(n,xc,w,wn):
In:=0:
for i from 1-ind to N do:
 In:=In+evalf(w[i]*subs(x=xc[i],f)):
end:
for i from 1 to N do:
 In:=In+evalf(w[i]*subs(x=-xc[i],f)):
end:

end proc:


ghnw:=proc(n,xc,w,wn)
local N,xp,epsil,ind,epso,delta,xn,u0,up1,i,err,t,um1,um2,suma,sumad,k,errod,up2,ntd,ws,hs,sumaz,nt,deri,ii,kk,expo;
description "Gauss-Hermite quadrature";
epsil:=evalf(10^(-Digits/4)):
N:=floor(n/2):
ind:=n-2*N:
epso:=evalf(10^(-Digits)):
xp:=0:
delta:=evalf((1+ind)/2*Pi/sqrt(2*n+1)):
xn:=evalf(xp+delta):
if ind>0 then
   u0:=0;up1:=1;
else
   u0:=1;up1:=0;
end if:
i:=0:
while i<N do
    i:=i+1:
    err:=1+epsil:
    while err>epsil do
        t:=delta;
        um1:=0;um2:=0;
        suma:=evalf(u0+up1*t):
        sumad:=up1;
        k:=-1:
        errod:=1;
        while errod>epso and k<Digits+20 do
             k:=k+1;
             up2:=evalf((xp^2-2*n-1)*u0+2*k*xp*um1+k*(k-1)*um2);
             um2:=um1;um1:=u0;u0:=up1;up1:=up2;
             ntd:=up2*t;
             sumad:=evalf(sumad+ntd);
             t:=t*delta/(k+2);
             nt:=up2*t;
             suma:=evalf(suma+nt);
             if i>1 and k>20 then
                 errod:=abs(evalf(nt/suma)):
             end:
       end do:
       u0:=suma;up1:=sumad;
       hs:=u0/up1:
       ws:=evalf(sqrt(2*n+1-xn^2));
       delta:=evalf(-arctan(ws*hs)/ws);
       err:=evalf(abs(delta)):
       xp:=xn;
       xn:=evalf(xp+delta);       
    end do:
    xc[i]:=xn:
    t:=delta;
    um1:=0;um2:=0;
    suma:=evalf(u0+up1*t):
    sumad:=up1;
    k:=-1;
    errod:=1;
    while errod>epso and k<Digits+20 do
          k:=k+1;
          up2:=evalf((xp^2-2*n-1)*u0+2*k*xp*um1+k*(k-1)*um2);
          um2:=um1;um1:=u0;u0:=up1;up1:=up2;
          ntd:=up2*t;
          sumad:=evalf(sumad+ntd);
          t:=t*delta/(k+2);
          nt:=up2*t;
          suma:=evalf(suma+nt);
          if i>1 and k>20  then
              errod:=abs(evalf(nt/suma));
          end;
    end do:
    u0:=suma;up1:=sumad;
    deri[i]:=sumad;
    delta:=evalf(Pi/sqrt(2*n+1-xn^2));
    xp:=xn;
    xn:=evalf(xp+delta);
end do:
sumaz:=0:
for i from 1 to N do;
    sumaz:=sumaz+2*xc[i]^2*exp(-xc[i]^2)/deri[i]^2;
end:
for i from 1 to N do;
    wn[i]:=evalf(sqrt(Pi)/2/sumaz/deri[i]^2);
    w[i]:=exp(-xc[i]^2)*wn[i]:
end:
ii:=1:
if ind>0 then
   ii:=0:
   xc[0]:=0:
   wn[0]:=evalf(sqrt(Pi)/2/sumaz):
   w[0]:=wn[0]:
end:

end proc;

