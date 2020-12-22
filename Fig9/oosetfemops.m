function p=oosetfemops(p)
% for triangular SKT cross-diffusion
% hence no K, it needs to be built in each step)
    gr=p.pdeo.grid;
    [~,M,~]=p.pdeo.fem.assema(gr,0,1,1); % assemble 'scalar' M
    p.mat.M=kron([[1,0];[0,1]],M); %build 2-component system M
    %% differentiation and convection matrices for the Jacobian
    p.mat.Dx=makeDx(p);
    p.mat.Kx=p.pdeo.fem.convection(gr,1);
end