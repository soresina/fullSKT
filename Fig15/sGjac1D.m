function Gu=sGjac1D(p,u) 
% jacobian for triSKT
u1=u(1:p.np);
u2=u(p.np+1,2*p.np);
par=u(p.nu+1:end);
d=par(1); d12=par(2); d21=par(3);
r1=par(4); r2=par(5); a1=par(6); a2=par(7); b1=par(8); b2=par(9);
s1=par(10); s2=par(11);

f1=(r1-a1*u1-b1*u2).*u1;
f2=(r2-b2*u1-a2*u2).*u2;

f1_u1=r1-2*a1*u1-b1*u2;
f1_u2=-b1*u1;
f2_u1=-b1*u2;
f2_u2=r2-b1*u1-2*a2*u2;

n=p.np;
Fu=[[spdiags(f1_u1,0,n,n),spdiags(f1_u2,0,n,n)];
    [spdiags(f2_u1,0,n,n),spdiags(f2_u2,0,n,n)]];

gr=p.pdeo.grid; 
u1t=gr.point2Center(u1); % interpolate to element centers
u2t=gr.point2Center(u2);

c=d+d12*u2t+s1*u1t;
cc=d12*u1t;	

c_u2=d12;
cc_u1=d12;

c2=d+d21*u1t;
cc2=d21*u2t;

c2_u1=d21;
cc2_u2=d21;

u1x=p.mat.Dx*u1; % 1st derivative as coefficient
u2x=p.mat.Dx*u2;
K1c=p.mat.Kx*spdiags(c_u2.*u2x,0,n,n); % first order derivative 
K1cc=p.mat.Kx*spdiags(cc_u1.*u1x,0,n,n); % first order derivative

K1c2=p.mat.Kx*spdiags(c2_u1.*u1x,0,n,n); % first order derivative 
K1cc2=p.mat.Kx*spdiags(cc2_u2.*u2x,0,n,n); % first order derivative 


[Kc,M1,~]=p.pdeo.fem.assema(gr,c,0,0); % K(u) for div(h(u)*grad v)  part  
[Kcc,~,~]=p.pdeo.fem.assema(gr,cc,0,0); % K(u) for div(h(u)*grad v)  part  
[Kc2,M2,~]=p.pdeo.fem.assema(gr,c2,0,0);
[Kcc2,~,~]=p.pdeo.fem.assema(gr,cc2,0,0); 

N=sparse(p.pdeo.grid.nPoints,p.pdeo.grid.nPoints);
KK=[[Kc Kcc];[Kcc2 Kc2]];
KK1=[[K1cc K1c];[K1cc2 K1c2]];
M=[[M1 N];[N M2]];
Gu=KK-KK1-M*Fu;  % putting it all together
end