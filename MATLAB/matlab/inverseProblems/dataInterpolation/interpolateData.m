clc
clear all
addpath(genpath('/scratch/waquino/code/Salinas/tools/matlab'));

%Arch1 side surface
dataSurface='side';
 nodesetID = 1;
 inputFile='input.mat';
 interpolateFile='delaminationRecon.mat';
 meshFileToTestInterpolation='delaminationRecon.exo';
 testInterpolation=1;

%%***********************************
load(interpolateFile);
[Fu, Fv, Fw] = constructInterpolant(inputFile);

nodeIndex = find(nsids == nodesetID);
 
targetNodes = eval(strcat('nsnod0',int2str(nodeIndex)));
%indexesFromSideNodes = find(node_num_map == sideNodes);
xside = x0(targetNodes);
yside = y0(targetNodes);
zside = z0(targetNodes);

dispX = Fu(xside, yside, zside);
dispY = Fv(xside, yside, zside);
dispZ = Fw(xside, yside, zside);

%Output files to disk
writeIputFilesForIP(targetNodes,node_num_map, dispX, dispY, dispZ);


%Add displacements to exodus file
%This code is just for testing interpolation
if (testInterpolation)
    db=exo_get(meshFileToTestInterpolation);
    db.NodalVars(1).Data(:)=0;
    db.NodalVars(2).Data(:)=0;
    db.NodalVars(3).Data(:)=0;
    for i=1:length(targetNodes)
        index=targetNodes(i);
        db.NodalVars(1).Data(index)=dispX(i);
        db.NodalVars(2).Data(index)=dispY(i);      
        db.NodalVars(3).Data(index)=dispZ(i);
    end;
    exo_put('interpolatedData-out.exo',db);
end
