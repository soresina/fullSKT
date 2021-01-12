function r=sG(p,u) 
% compute pde-part of residual
u1=u(1:p.np); % extract the first component
u2=u(p.np+1:2*p.np); % extract the second component
par=u(p.nu+1:end); % extract parameters
d=par(1); d12=par(2); d21=par(3);
r1=par(4); r2=par(5); a1=par(6); a2=par(7); b1=par(8); b2=par(9);
s1=par(10); s2=par(11);

f1=(r1-a1*u1-b1*u2).*u1;
f2=(r2-b2*u1-a2*u2).*u2;

gr=p.pdeo.grid;
%ut=(p.mat.p2c*u)'; % interpolate to element centers
u1t=gr.point2Center(u1);
u2t=gr.point2Center(u2);

c=d+d12*u2t+2*s1*u1t;
cc=d12*u1t;

c2=d+d21*u1t+2*s2*u2t;
cc2=d21*u2t;

[Kc,~,F1]=p.pdeo.fem.assema(gr,c,0,f1);
[Kcc,~,~]=p.pdeo.fem.assema(gr,cc,0,0);
[Kc2,~,F2]=p.pdeo.fem.assema(gr,c2,0,f2);
[Kcc2,~,~]=p.pdeo.fem.assema(gr,cc2,0,0);
%N=sparse(p.pdeo.grid.nPoints,p.pdeo.grid.nPoints);
p.mat.K=[[Kc Kcc];[Kcc2 Kc2]];
F=[F1;F2];
r=p.mat.K*[u1;u2]-F;
end