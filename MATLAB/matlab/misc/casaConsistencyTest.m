%CASACONSISTENCYTEST 3. check consistency of matrices and vectors
% addpath([pwd,'/handoff']); % only once uncomment
clear
p = 8;
% name = 'getphi_gid_';
% [localMap,globalRow]=parallelmap(name,p);
load('casamaps.mat');
% $ tail -n 2 Kssr*.m
suborder=[36432,50196,36348,36348,43554,37578,38709,38034];
for i = 1:p,
  assert(suborder(i)==size(globalRow{i},2));
end
order = getLastGlobalRow(globalRow);
load('rbm'); % casaDisp
numdof = size(rbm,1);
assert( numdof == order );
disp('passed');