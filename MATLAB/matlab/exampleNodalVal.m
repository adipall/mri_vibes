addpath /scratch/dmday/code/Salinas/tools/matlab/loadExodusIntoMatlab
name = 'preload.exo';
skipNodalVars = 6;
fexo=exo_rd(name,skipNodalVars);
clear name skipNodalVars
c = fexo.Nodes.Coordinates;
n = size(c,1);
temperature =  ones(n,1)*927- c(:,1)*42.45; 
fexo.NodalVars(1).Name = 'temperature';
fexo.NodalVars(1).Data = temperature;
filename = 'preload_thermal.exo';
exo_put(filename,fexo);
