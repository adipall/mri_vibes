%% this script takes transient history file output (synthetic data) and imprints it onto 
%% the original exodus file as nodeset data. The acoustic and structural nodesets
%% are specified below


% this is a history file that has the synthetic data
load hbeam_transient-out.mat

dataTimes = time;

% these are the nodeset numbers of the acoustic and structural nodesets where history data was written
acousticNodeSetID = 102;
structuralNodeSetID = 25;

% this is the original (large) mesh file
baseExodusFile = 'hbeam.exo';
new=exo_get(baseExodusFile)

nsetids = [new.Nodesets(:).ID];

% these indicies come from the large (original) exodus file. They are used only to index into the 'new' struct
acousticNodeSetIndex = 1;
structuralNodeSetIndex = 2;

%these arrays come from the .frq or history file that contains the data. They are used only to extract
% pressure and displacement data from the history or frq arrays
structuralNodes = nsnod01;
acousticNodes = nsnod02;

dispxdata=nvar01(structuralNodes,:);
dispydata=nvar02(structuralNodes,:);
dispzdata=nvar03(structuralNodes,:);
pressuredata=nvar07(acousticNodes,:);

new.NodesetVars(structuralNodeSetIndex).NodesetVarData(1).Name = 'Dispx';
new.NodesetVars(structuralNodeSetIndex).NodesetVarData(1).Data = dispxdata;

new.NodesetVars(structuralNodeSetIndex).NodesetVarData(2).Name = 'Dispy';
new.NodesetVars(structuralNodeSetIndex).NodesetVarData(2).Data = dispydata;

new.NodesetVars(structuralNodeSetIndex).NodesetVarData(3).Name = 'Dispz';
new.NodesetVars(structuralNodeSetIndex).NodesetVarData(3).Data = dispzdata;

new.NodesetVars(structuralNodeSetIndex).NodesetVarData(4).Name = 'APressure';
new.NodesetVars(structuralNodeSetIndex).NodesetVarData(4).Data = [];

new.NodesetVars(acousticNodeSetIndex).NodesetVarData(1).Name = 'Dispx';
new.NodesetVars(acousticNodeSetIndex).NodesetVarData(2).Name = 'Dispy';
new.NodesetVars(acousticNodeSetIndex).NodesetVarData(3).Name = 'Dispz';
new.NodesetVars(acousticNodeSetIndex).NodesetVarData(1).Data = [];
new.NodesetVars(acousticNodeSetIndex).NodesetVarData(2).Data = [];
new.NodesetVars(acousticNodeSetIndex).NodesetVarData(3).Data = [];

% we need to convert pressure to velocity potential
% note that this assumes that the pressure data is padded with some zeros at the beginning
numpoints = size(acousticNodes,1);
for inode = 1:numpoints
   velpotdata(inode,:) = cumtrapz(dataTimes,pressuredata(inode,:));
end

new.NodesetVars(acousticNodeSetIndex).NodesetVarData(4).Name = 'APressure';
new.NodesetVars(acousticNodeSetIndex).NodesetVarData(4).Data = velpotdata;

new.Time = dataTimes;

exo_put('hbeam_transient_new.exo',new)

