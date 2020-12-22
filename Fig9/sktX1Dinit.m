function p=sktX1Dinit(p,lx,nx,par)
% (generic) init routine for the
% skt cross-diffusion system
% u_t=Lap((d+d12 v) u)+(r1-a1 u -b1 v)u
% v_t=Lap((d+d21 u) v)+(r2-b2 u -a1 v)v
    %% setting standard parameters and screenlayout
    p=stanparam(p); % infuses p with standard parameter settings
    screenlayout(p); % open, clear and arrange the common figures
    %% special parameters related to this model
    % basics
    p.nc.neq=2; % number of equations in the model
    p.nc.nsteps=200;
    p.sw.sfem=-1; % type of numerical calculation, here OOPDE
    p.sw.bifcheck=2;
    p.sw.foldcheck=1;
    p.sw.jac=0; % use numerical jacobian
    p.sw.spjac=0; % use analytical Jacobian for spectral point cont (fold cont)
    p.plot.pcmp=1; % plotsol plots the 2nd component
    % description of the model
    p.fuha.sG=@sG; % the model itself
    p.fuha.sGjac=@sGjac1D; % the Jacobian of the model  
    p.fuha.outfu=@sktX1Dbra; % output quantities
    %% domain and mesh
    p.pdeo=stanpdeo1D(lx,2*lx/nx); % mesh [-lx ,lx], max mesh pt 2* lx/r
    bc=p.pdeo.grid.neumannBC('0');
    p.pdeo.grid.makeBoundaryMatrix(bc);
    p.np=p.pdeo.grid.nPoints ; % number of meshpoints
    p.nu=p.np*p.nc.neq; % number of unknowns (=2*( mesh points ), as 3 components )
    p=setfemops(p); % compute FEM - operators
    %% bifurcation parameter , continuation basics and first guess for solution
    p.nc.ilam=1; % primary bifurcation parameter located at p.u(p.np+p.nc.ilam )
    p.nc.lammin=1e-5; % lower bound for primary parameter during continuation
    p.nc.lammax=0.41; % upper bound for primary parameter during continuation
    p.sol.xi=1/p.nu; % weight in arclength - continuation
    p.sol.ds=-0.01; % starting stepsize
    p.nc.dsmax=1e-3; % maximal stepsize
    p.nc.dsmin=1e-6; % minimal stepsize
    %p.nc.ngen=100;
    % construction the trivial solution
    % initial guess
    r1=par(4); r2=par(5); a1=par(6); a2=par(7); b1=par(8); b2=par(9);
    v=(a1*r2-r1*b2)/(a1*a2-b1*b2);
    u=(r1-b1*v)/a1; 
    p.u=[u*ones(p.np,1);v*ones(p.np,1);par]; % initial solution guess with parameters
    %% Plot
    p.plot.auxdict ={'d','d12','d21','r1','r2','a1','a2','b1','b2'}; % parameter names 

end