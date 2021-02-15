%function aCell = loadConstraintMatrix(p)
% See also loadMatrix.m
% example casaConstraint.m
function aCell = loadConstraintMatrix(p)
prefix = 'matrix=ConstraintMatrix_';
aCell = cell(p,1);
suffix = '();';
matrix=[];
for i = 1:p,
    eval([prefix,int2str(i-1),suffix]);
    aCell{i}=matrix;
end