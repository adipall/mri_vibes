% addpath /scratch/dmday/code/Salinas/tools/matlab/loadExodusIntoMatlab
name = 'nail.exo';
fexo=exo_rd(name);
nvars=length(fexo.ElemVars(1,:));
eine = ones(6,1);
zilch = zeros(6,1);
scale = -5.e+9;
for i=1:27,
    gp = 6*(i-1);
    fexo.ElemVars(1,gp+1).Data = zilch;
    fexo.ElemVars(1,gp+2).Data = eine*scale;
    fexo.ElemVars(1,gp+3).Data = zilch;
    fexo.ElemVars(1,gp+4).Data = zilch;
    fexo.ElemVars(1,gp+5).Data = zilch;
    fexo.ElemVars(1,gp+6).Data = zilch;
end
exo_put('preload.exo',fexo);
