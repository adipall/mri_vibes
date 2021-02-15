%function aCell = loadMatrix(p)
% See also loadMaps
% example
function aCell = loadMatrix(p)
prefix = 'matrix=Kssr_';
aCell = cell(p,1);
suffix = '();';
matrix=[];
for i = 1:p,
    eval([prefix,int2str(i-1),suffix]);
    matrix = matrix + tril(matrix,-1)';
    aCell{i}=matrix;
end
