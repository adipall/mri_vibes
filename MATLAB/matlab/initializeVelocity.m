% Sierra SD spatially dependent fields 
% Next 3 lines go in Sierra SD text input file
%%INITIAL-CONDITIONS
%%   velocity from_file //displacement vanishes
%%END
%
% This Matlab script adds a nodal variable to a Genesis file.
% It has not been tested with an Exodus file that already
% contains nodal variables (as far as I know).
% 
clear
% Here the name of the output file is ...
filename = 'initialized.exo';
% If re-running this, then the previous file must be removed
% first.  This is done from the shell, not from Matlab.
% $ rm initialized.exo
% The name of Matt's Genesis file is ...
name='solid_3d_simply_supported.g';
%-------- only set velocity below should change -------
skipNodalVars=6;
fexo=exo_rd(name,skipNodalVars);
clear name skipNodalVars
coords=fexo.Nodes.Coordinates;
n=size(coords,1);
velocity=zeros(n,6);
velocity(:,2)=sin(pi*coords(:,1)).*sin(pi*coords(:,3));%set velocity
clear coords n
node_names={'VelX','VelY','VelZ','VelRX','VelRY','VelRZ'};
time=0;
for i=1:6,
    fexo = AddNodalVar(fexo,node_names{i},velocity(:,i),time);
end
exo_put(filename,fexo);
