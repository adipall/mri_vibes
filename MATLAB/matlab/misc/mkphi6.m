function phi=mkphi6(nsteps,nvar01,nvar02,nvar03,nvar04,nvar05,nvar06)
%function phi=mkphi6(nsteps,nvar01,nvar02,nvar03,nvar04,nvar05,nvar06)
%compute the phi matrix from nvar01 through nvar03 in the
% example
%   load('foo.mat');  % created from exo2mat foo.exo
%   phi = mkphi6(nsteps,nvar01,nvar02,nvar03,nvar04,nvar05,nvar06);
%   phi is defined on the structural set, with 6 dof per node.
nnodes = size(nvar01,1);
phi=zeros(nnodes*6,nsteps);
tmp=(0:nnodes-1)*6;
phi(tmp+1,:)=nvar01;
phi(tmp+2,:)=nvar02;
phi(tmp+3,:)=nvar03;
phi(tmp+4,:)=nvar04;
phi(tmp+5,:)=nvar05;
phi(tmp+6,:)=nvar06;
