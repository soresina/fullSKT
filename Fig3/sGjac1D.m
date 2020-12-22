function Gu=sGjac1D(p,u) 
% jacobian for triSKT
u1=u(1:p.np);
u2=u(p.np+1,2*p.np);
par=u(p.nu+1:end);
d=par(1); alpha=par(2);
r1=par(3); r2=par(4); a1=par(5); a2=par(6); b1=par(7); b2=par(8);

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

c=d+alpha*u2t;
cc=alpha*u1t;	

c_u2=alpha;
cc_u1=alpha;

u1x=p.mat.Dx*u1; % 1st derivative as coefficient
u2x=p.mat.Dx*u2;
K1c=p.mat.Kx*spdiags(c_u2.*u2x,0,n,n); % first order derivative 
K1cc=p.mat.Kx*spdiags(cc_u1.*u1x,0,n,n); % first order derivative 

[Kc,M1,~]=p.pdeo.fem.assema(gr,c,0,f1); % K(u) for div(h(u)*grad v)  part  
[Kcc,~,~]=p.pdeo.fem.assema(gr,cc,0,f1); % K(u) for div(h(u)*grad v)  part  
[K,M2,~]=p.pdeo.fem.assema(gr,1,0,f2);
N=sparse(p.pdeo.grid.nPoints,p.pdeo.grid.nPoints);
KK=[[Kc Kcc];[N d*K]];
KK1=[[K1cc K1c];[N N]];
M=[[M1 N];[N M2]];
Gu=KK-KK1-M*Fu;  % putting it all together
end