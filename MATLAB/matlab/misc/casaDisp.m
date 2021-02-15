% 2. convert mode shapes from structural set to analysis set
% addpath([pwd,'/handoff']); % only once uncomment
clear
p = 8;
% name = 'getcorrectedphi_gid_';
% [localMap,globalRow]=parallelmap(name,p);
load('casamaps.mat');
order = getLastGlobalRow(globalRow);
% use SierraSD to compute the rigid body modes
% $ epu -auto 1/casa-rigid.par.8.0
% $ mv casa-rigid.par casagrbm.exo
% $ exo2mat casagrbm.exo
load('casagrbm.mat');% nvar01,nvar02,nvar03
dofpernode = 3;
phi=mkphi3(nvar01,nvar02,nvar03);
numdof = size(phi,1);
assert( numdof >= order );
numModes = size(phi,2);
rbm = zeros(order,numModes);

for i = 1:p,
    row = globalRow{i}+1;
    ordinal = dofpernode*(localMap{i}(:,2)-1)+localMap{i}(:,3);
    rbm(row,:) = phi(ordinal,:);
end
assert( size(rbm,1) == order );
save('rbm.mat','rbm');
