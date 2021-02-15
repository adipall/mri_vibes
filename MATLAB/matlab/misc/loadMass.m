%function aCell = loadMass(p)
% See also loadMatrix
% example loadMass
function aCell = loadMass(p)
prefix = 'matrix=Mssr_';
aCell = cell(p,1);
suffix = '();';
matrix=[];
for i = 1:p,
    eval([prefix,int2str(i-1),suffix]);
    matrix = matrix + tril(matrix,-1)';
    aCell{i}=matrix;
end