%CASAMASSTOO optional, e.g. check R'MR=eye
addpath([pwd,'/handoff']); % only once uncomment
clear
p = 8;
submatrix = loadMass(p);
load('casamaps.mat');
order=getLastGlobalRow(globalRow);
mass = sparse(order,order);
for i = 1:p,
  row = globalRow{i}+1;
  mass(row,row) = mass(row,row) + submatrix{i};
end
clear submatrix row maxrow
