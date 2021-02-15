%function C = loadGids(root,p)
% See also dofmap
% example
% overlaps.m
function C = loadGids(root,p)
C = cell(p,1);
for i = 1:p,
    eval([root,int2str(i-1)]),
    C{i} = gid;
end
