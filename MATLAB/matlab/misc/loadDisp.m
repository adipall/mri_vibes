%function aCell = loadDisp(p)
% See also loadMaps
% example
function aCell = loadDisp(p)
disp('Warning: disp files are indexed in a way that makes no sense to developers');
prefix = 'phi=getphi_Disp';
aCell = cell(6*p,1);
suffix = '();';
phi=[];
for mode=1:6,
  for i = 1:p,
    eval([prefix,int2str(mode),'_',int2str(i-1),suffix]);
    j = 6*(mode-1)+i;
    aCell{j}=phi;
  end
end
