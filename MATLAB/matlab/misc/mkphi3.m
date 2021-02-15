% assemble the phi matrix from nvar01 through nvar03
% Hypothesis. entirely solid model: nvar04 = nvar05 = nvar06 = 0
% example   
%   load('casagrbm.mat');  % created from exo2mat foo.exo
%   rbm = mkphi3(nvar01,nvar02,nvar03);
%   save('rbm.mat','rbm');
% function phi=mkphi3(nvar01,nvar02,nvar03)
function phi=mkphi3(nvar01,nvar02,nvar03)
numNodes = size(nvar01,1);  % N
numModes = size(nvar01,2);  % M
phi=zeros(3*numNodes,numModes);
phi(1:3:end,:)=nvar01;
phi(2:3:end,:)=nvar02;
phi(3:3:end,:)=nvar03;
