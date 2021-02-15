% example loadGids
% See also loadGids.m
% addpath([pwd,'/example']);
root = 'tet4_small_statics_gid_';
p = 6;
C = loadGids(root,p);
pm1 = p-1;
for i = 1:pm1,
    for j=i+1:p,
        cap = intersect(C{i},C{j});
        if isempty(cap)
            disjoint = [i,j];
        end
    end
end