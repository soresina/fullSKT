function p=sktX1Dinit(p,lx,nx,par)
% (generic) init routine for the
% triangular skt cross-diffusion system
% u_t=Lap((d+alpha v) u)+(r1-a1 u -b1 v)u
% v_t=d Lap v           +(r2-b2 u -a1 v)v
    %% setting standard parameters and screenlayout
    p=stanparam(p); % infuses p with standard parameter settings
    screenlayout(p); % open, clear and arrange the common figures
    %% special parameters related to this model
    % basics
    p.nc.neq=2; % number of equations in the model
    p.nc.nsteps=100;
    p.sw.foldcheck=1;
    p.sw.bifcheck=2;
    p.sw.sfem=-1; % type of numerical calculation, here OOPDE
    p.sw.jac=0; % use numerical jacobian
    p.sw.spjac=0; % use analytical Jacobian for spectral point cont (fold cont)
    p.plot.pcmp=2; % plotsol plots the 2nd component
    % description of the model
    p.fuha.sG=@sG; % the model itself
    %p.fuha.sGjac=@sGjac1D; % the Jacobian of the model  
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
    p.nc.lammin=1e-4; % lower bound for primary parameter during continuation
    p.nc.lammax=0.04; % upper bound for primary parameter during continuation
    p.sol.xi=1/p.nu; % weight in arclength - continuation
    p.sol.ds=-0.0001; % starting stepsize
    p.nc.dsmax=5e-4; % maximal stepsize
    p.nc.dsmin=1e-7; % minimal stepsize
    % construction the trivial solution
    % initial guess
    r1=par(3); r2=par(4); a1=par(5); a2=par(6); b1=par(7); b2=par(8);
    v=(a1*r2-r1*b2)/(a1*a2-b1*b2);
    u=(r1-b1*v)/a1; 
    p.u=[u*ones(p.np,1);v*ones(p.np,1);par]; % initial solution guess with parameters
    %% Plot
    p.plot.auxdict ={'d','d12','r1','r2','a1','a2','b1','b2'}; % parameter names 

end