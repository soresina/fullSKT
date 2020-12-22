function r=sG(p,u) 
% compute pde-part of residual
u1=u(1:p.np); % extract the first component
u2=u(p.np+1:2*p.np); % extract the second component
par=u(p.nu+1:end); % extract parameters
d=par(1); alpha=par(2);
r1=par(3); r2=par(4); a1=par(5); a2=par(6); b1=par(7); b2=par(8);

f1=(r1-a1*u1-b1*u2).*u1;
f2=(r2-b2*u1-a2*u2).*u2;

gr=p.pdeo.grid;
%ut=(p.mat.p2c*u)'; % interpolate to element centers
u1t=gr.point2Center(u1);
u2t=gr.point2Center(u2);

c=d+alpha*u2t;
cc=alpha*u1t;
[Kc,~,F1]=p.pdeo.fem.assema(gr,c,0,f1);
[Kcc,~,~]=p.pdeo.fem.assema(gr,cc,0,f1);
[K,~,F2]=p.pdeo.fem.assema(gr,1,0,f2);
N=sparse(p.pdeo.grid.nPoints,p.pdeo.grid.nPoints);
p.mat.K=[[Kc Kcc];[N d*K]];
F=[F1;F2];
r=p.mat.K*[u1;u2]-F;
end