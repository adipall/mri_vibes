%% this script takes directfrf history file output (synthetic data) and imprints it onto 
%% the original exodus file as nodeset data. The acoustic and structural nodesets
%% are specified below

% this is a history file that has the synthetic data
load hbeam-out.mat

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
imagdispxdata=nvar02(structuralNodes,:);
dispydata=nvar03(structuralNodes,:);
imagdispydata=nvar04(structuralNodes,:);
dispzdata=nvar05(structuralNodes,:);
imagdispzdata=nvar06(structuralNodes,:);
pressuredata=nvar13(acousticNodes,:);
imagpressuredata=nvar14(acousticNodes,:);

new.NodesetVars(structuralNodeSetIndex).NodesetVarData(1).Name = 'Dispx';
new.NodesetVars(structuralNodeSetIndex).NodesetVarData(1).Data = dispxdata;
new.NodesetVars(structuralNodeSetIndex).NodesetVarData(2).Name = 'ImagDispx';
new.NodesetVars(structuralNodeSetIndex).NodesetVarData(2).Data = imagdispxdata;

new.NodesetVars(structuralNodeSetIndex).NodesetVarData(3).Name = 'Dispy';
new.NodesetVars(structuralNodeSetIndex).NodesetVarData(3).Data = dispydata;
new.NodesetVars(structuralNodeSetIndex).NodesetVarData(4).Name = 'ImagDispy';
new.NodesetVars(structuralNodeSetIndex).NodesetVarData(4).Data = imagdispydata;

new.NodesetVars(structuralNodeSetIndex).NodesetVarData(5).Name = 'Dispz';
new.NodesetVars(structuralNodeSetIndex).NodesetVarData(5).Data = dispzdata;
new.NodesetVars(structuralNodeSetIndex).NodesetVarData(6).Name = 'ImagDispz';
new.NodesetVars(structuralNodeSetIndex).NodesetVarData(6).Data = imagdispzdata;

new.NodesetVars(structuralNodeSetIndex).NodesetVarData(7).Name = 'APressure';
new.NodesetVars(structuralNodeSetIndex).NodesetVarData(7).Data = [];
new.NodesetVars(structuralNodeSetIndex).NodesetVarData(8).Name = 'ImagAPressure';
new.NodesetVars(structuralNodeSetIndex).NodesetVarData(8).Data = [];

new.NodesetVars(acousticNodeSetIndex).NodesetVarData(1).Name = 'Dispx';
new.NodesetVars(acousticNodeSetIndex).NodesetVarData(2).Name = 'ImagDispx';
new.NodesetVars(acousticNodeSetIndex).NodesetVarData(3).Name = 'Dispy';
new.NodesetVars(acousticNodeSetIndex).NodesetVarData(4).Name = 'ImagDispy';
new.NodesetVars(acousticNodeSetIndex).NodesetVarData(5).Name = 'Dispz';
new.NodesetVars(acousticNodeSetIndex).NodesetVarData(6).Name = 'ImagDispz';
new.NodesetVars(acousticNodeSetIndex).NodesetVarData(1).Data = [];
new.NodesetVars(acousticNodeSetIndex).NodesetVarData(2).Data = [];
new.NodesetVars(acousticNodeSetIndex).NodesetVarData(3).Data = [];
new.NodesetVars(acousticNodeSetIndex).NodesetVarData(4).Data = [];
new.NodesetVars(acousticNodeSetIndex).NodesetVarData(5).Data = [];
new.NodesetVars(acousticNodeSetIndex).NodesetVarData(6).Data = [];

% note we need to convert pressure to velocity potential
temprealpressure = pressuredata;
tempimagpressure = imagpressuredata;
for ifreq=1:length(dataTimes)
   freq = dataTimes(ifreq);
   velpotdata(:,ifreq) = tempimagpressure(:,ifreq)/(2*pi*freq);
   imagvelpotdata(:,ifreq) = -temprealpressure(:,ifreq)/(2*pi*freq);
end

new.NodesetVars(acousticNodeSetIndex).NodesetVarData(7).Name = 'APressure';
new.NodesetVars(acousticNodeSetIndex).NodesetVarData(7).Data = velpotdata;
new.NodesetVars(acousticNodeSetIndex).NodesetVarData(8).Name = 'ImagAPressure';
new.NodesetVars(acousticNodeSetIndex).NodesetVarData(8).Data = imagvelpotdata;

new.Time = dataTimes;

exo_put('hbeam_new.exo',new)

