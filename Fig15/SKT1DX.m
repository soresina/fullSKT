% SKT - triangular cross diffusion system 
% Lap[(d1+d12 v+d11 u)u]+(r1-a1u-b1v)u=0
% Lap[(d1+d21 u+d22 v)v]+(r2-b2u-a2u)v=0
% !! d1=d2=d bifurcation parameter
clear
clc
close all
%% 1 - creating a clear working space and setting pde2path
cFolPath='/Users/cinziasoresina/Google Drive/UniversityAndWorks/TUM/pde2path/pde2path/myworks/SKT_GitHub/Fig15';
cd('/Users/cinziasoresina/Google Drive/UniversityAndWorks/TUM/pde2path/pde2path')

setpde2path
edtSvc=com.mathworks.mlservices.MLEditorServices; 
edtList=edtSvc.getEditorApplication.getOpenEditors.toArray;
for k=1:length(edtList)
    [~, fname]=fileparts(char(edtList(k).getLongName.toString));
    edt.(fname)=edtList(k);
end
edt.setpde2path.closeNoPrompt % close setpde2path.m
cd(cFolPath)
close all
keep pphome cFolPath
clc
%% 2 - initialising the problem
p=[];
% 'd','d12','d21','r1','r2','a1','a2','b1','b2','d11','d22'
par=[0.04,3,0,5,2,3,3,1,1,0,0.03]';

lx=0.5;  % domain (interval)
p=sktX1Dinit(p,lx,50,par); 
p.Om=2*lx;
p=setfn(p,'hom_branch'); % current directory
%% 3 - contiunation of the trivial branch
p=cont(p,100); % continuation for a maximum of x steps
%% 4 - switch to periodic branch and continuation with fold detection
p=swibra('hom_branch','bpt1','bpt1',1e-2); p.nc.dsmax=1e-2; p=cont(p,100);
%%
nmax=200;
p=swibra('hom_branch','bpt2','bpt2_u',1e-2);  p.nc.dsmax=1e-2; p=cont(p,nmax);
p=swibra('hom_branch','bpt2','bpt2_d',-1e-2);  p.nc.dsmax=1e-2; p=cont(p,nmax);
%%
p=swibra('hom_branch','bpt3','bpt3_u',1e-2); p.nc.dsmax=1e-2; p=cont(p,nmax);
p=swibra('hom_branch','bpt3','bpt3_d',-1e-2);  p.nc.dsmax=1e-2; p=cont(p,nmax); 
%%
p=swibra('bpt2_d','bpt1','b2d_u',1e-2); p=cont(p,nmax);
%%
p=swibra('bpt3_d','bpt1','b3d_u',1e-2); p.nc.lammin=-10; p=cont(p,nmax);
%% Postprocessing, first plot BifDiagram 
nfig=151;
figure(nfig);
clf(nfig); 
cmp=0;
box on
hold on
plotbra('hom_branch',nfig,cmp,'cl','k'); 
plotbra('bpt1',nfig,cmp,'cl','b');
plotbra('bpt2_d',nfig,cmp,'cl','r');
plotbra('bpt2_u',nfig,cmp,'cl','r');
plotbra('bpt3_d',nfig,cmp,'cl','g'); 
plotbra('bpt3_u',nfig,cmp,'cl','g');
plotbra('b2d_u',nfig,cmp,'cl','m');
plotbra('b3d_u',nfig,cmp,'cl','m');

xlabel('d')
ylabel('')
%axis([0 0.04 0 0.5])
axis([0 0.04 1.5 1.65])
