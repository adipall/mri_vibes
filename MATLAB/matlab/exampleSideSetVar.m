% addpath /scratch/dmday/code/Salinas/tools/matlab/loadExodusIntoMatlab
name = 'frustum.exo';
fexo=exo_rd(name);
varname = 'pressure';
side = 2;
numberSides = size(fexo.Sidesets(1,side).Elements,1)
data = ones(numberSides,1)*1.e7;
time = 0.;
fexo = AddSidesetVar(fexo,varname,side,data,time);
filename = 'frustum2.exo';
exo_put(filename,fexo);
