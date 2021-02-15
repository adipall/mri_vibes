% 4. Compute force = K*rbm, which should vanish
% addpath([pwd,'/handoff']); % only once uncomment
clear
p = 8;
% name = 'getphi_gid_';
% [localMap,globalRow]=parallelmap(name,p);
load('casamaps.mat');
submatrix = loadMatrix(p);
order = getLastGlobalRow(globalRow);
a = sparse(order,order);
for i = 1:p,
  row = globalRow{i}+1;
  a(row,row) = a(row,row) + submatrix{i};
end
assert( size(a,1) == order );
clear i row submatrix
load('rbm.mat'); % casaDisp
force = a*rbm;

% norm K = 7e+7,  norm( rbm ) = 900,
% norm force = 4 e-4 is great.