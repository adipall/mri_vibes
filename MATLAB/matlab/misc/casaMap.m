% 1. Generate analysis and structural set maps
addpath([pwd,'/handoff']); % only once uncomment
clear
name = 'getcorrectedphi_gid_';
p = 8;
[localMap,globalRow]=parallelmap(name,p);
save('casamaps.mat','localMap','globalRow');