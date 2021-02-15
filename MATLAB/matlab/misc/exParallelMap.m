% example 
% addpath([pwd,'/example']); % only once uncomment
clear
name = 'tet4_small_statics_gid_';
p = 6;
% See also overlaps.m, loadGids.m
% addpath([pwd,'/example']);% only once uncomment
[localMap,globalRow]=parallelmap(name,p);
submatrix = loadMatrix(p);
order = getLastGlobalRow(globalRow);
a = sparse(order,order);
for i = 1:p,
  row = globalRow{i}+1;
  a(row,row) = a(row,row) + submatrix{i};
end
clear i name row submatrix
