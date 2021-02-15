%EXTRIPLET example triplet
% See also overlaps.m, loadGids.m
% addpath([pwd,'/example']);% only once uncomment
p = 6;
root = 'tet4_small_statics_gid_';
g = loadGids(root,p);
m = loadMaps(p);
aMap = triplet(g{1},m{1});
