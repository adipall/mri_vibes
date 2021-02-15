%% this script takes transient history file output (synthetic data) and imprints it onto 
%% the original exodus file as nodeset data. The acoustic and structural nodesets
%% are specified below


% this is a history file that has the synthetic data
load ellipseInSquare-out.mat

dataTimes = time;

% these are the nodeset numbers of the acoustic and structural nodesets where history data was written
acousticNodeSetID = 1;

% this is the original (large) mesh file
baseExodusFile = 'ellipseInSquare.e';
new=exo_get(baseExodusFile)

nsetids = [new.Nodesets(:).ID];

% these indicies come from the large (original) exodus file. They are used only to index into the 'new' struct
structuralNodeSetIndex = 2;
acousticNodeSetIndex = 1;

%these arrays come from the .frq or history file that contains the data. They are used only to extract
% pressure and displacement data from the history or frq arrays
acousticNodes = nsnod01;

pressuredata=nvar01(acousticNodes,:);

new.NodesetVars(structuralNodeSetIndex).NodesetVarData(1).Name = 'APressure';
new.NodesetVars(structuralNodeSetIndex).NodesetVarData(1).Data = [];

%integrate in time to get velocity potential
VelPotential = pressuredata*0; % just to get dimensions right
deltaT = time(2)-time(1);
for i=1:length(time)-1
    VelPotential(:,i+1) = VelPotential(:,i) + 0.5*deltaT*(pressuredata(:,i) + pressuredata(:,i+1));
end

new.NodesetVars(acousticNodeSetIndex).NodesetVarData(1).Name = 'APressure';
new.NodesetVars(acousticNodeSetIndex).NodesetVarData(1).Data = VelPotential;

new.Time = dataTimes;

exo_put('ellipseInSquare_data.e',new)

