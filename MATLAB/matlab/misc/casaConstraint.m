% 5.A. load constraint matrices
% addpath([pwd,'/handoff/geometric']); % only once uncomment
clear
p = 8;
submatrix = loadConstraintMatrix(p);
numrow = size(submatrix{1},1);
for i = 1:p,
  assert( numrow == size(submatrix{i},1) );
end
suborder=[36432,50196,36348,36348,43554,37578,38709,38034]';
load('casamaps.mat');
order = getLastGlobalRow(globalRow);
constraint= sparse(numrow,order);
for i = 1:p,
  numcol = suborder(i);
  assert(  size( submatrix{i},2)<= numcol );
  [rowC,col,val]=find(submatrix{i});
  rowK = globalRow{i}+1;
  assert( size(rowK,2) == numcol );
  %constraint(:,rowK)=constraint(:,rowK)+sparse(rowC,col,val,numrow,numcol);
  constraint(:,rowK)=sparse(rowC,col,val,numrow,numcol);
end
clear rowK rowC suborder numcol val i col
load('rbm.mat'); % rbm=mkphi3(...)
r = constraint*rbm;
% clear submatrix order
