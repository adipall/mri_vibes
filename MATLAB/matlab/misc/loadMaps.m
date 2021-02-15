%function aCell = loadMaps(p)
% See also dofmap
% example
function aCell = loadMaps(p)
prefix = 'aSet=FetiMap_a_';
aCell = cell(p,1);
suffix = '();';
aSet =[];
for i = 1:p,
    eval([prefix,int2str(i-1),suffix]),
    aCell{i} = aSet;
end