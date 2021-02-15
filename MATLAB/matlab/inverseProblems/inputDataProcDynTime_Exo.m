function inputDataProcDynTime_Exo(nodesetId)

if (nargin < 1)
   nodesetId=1;
end

clc
%clear all

format long
%****************************INPUT PARAMETERS******************

% the file below contains the results from a forward run
inputfile ='cube_nodeset-out.exo';
% Measured dofs
X = 1;
Y = 1;
Z = 1;
noise = 0.;
% Measured dof names
stringNames={'DispX','DispY','DispZ'};
% stringNames={'AccX','AccY','AccZ'};
%Outputfile names
truthtablefile = './ttable.txt';
datafile = './data.txt';
%These paths need to be set to use exo_rd -- lets you get variable by stringNames
addpath('/scratch/lbmunda/code/Salinas/tools/matlab/loadExodusIntoMatlab/') %for exo_rd
addpath('/scratch/lbmunda/code/Salinas/tools/matlab') %for FEMesh

%***************************************************************
flag=0;
exodb = exo_rd(inputfile,flag);

nodeset = 0;
for i=1:length(exodb.Nodesets)
    if exodb.Nodesets(i).ID == nodesetId
        nodeset = i;
    end
end

if nodeset == 0
    disp('Nodeset ID not Found');
end
    
nodes = exodb.Nodesets(nodeset).Nodes;
locations = length(nodes);

steps = length(exodb.Time);

if (steps > 1)
    dt = exodb.Time(2) - exodb.Time(1);
else
    dt = exodb.Time(1);
end


%Construct the data file
did = fopen(datafile,'w');
  
fprintf(did, '%d  %d %d\n', 3*locations, steps, dt);

node_num_map = exodb.Nodes.NodeNumMap;
num_nodes_in_mesh = length(exodb.Nodes.Coordinates);
if isempty(node_num_map)
   for i=1:num_nodes_in_mesh
      node_num_map(i) = i;
   end
end

for d=1:length(stringNames)
    stringIndex=exodb.NodalVars.isNodalVar(stringNames(d));
    if stringIndex == 0,
        disp('stringNames(d) is not a NodalVar');
    end
    assert(stringIndex ~= 0);
    for p=1:locations
        for s=1:steps
            displ = exodb.NodalVars(stringIndex).Data(nodes(p), s);
            if (noise > 0)
                displ = displ*(1 + noise*randn());
                exodb.NodalVars(stringIndex).Data(nodes(p), s) = displ;
            end
            fprintf(did, '%17.12e     ', displ);
        end
        fprintf(did,'\n');
    end
end


%Construct the truth table
tid = fopen(truthtablefile,'w');
fprintf(tid, '%d\n', locations);

for i=1:locations
    mapped_node = node_num_map(nodes(i));
    fprintf(tid,'%d  ', mapped_node);
    fprintf(tid, '%d %d %d\n', X, Y, Z);
end
fclose(tid);

fclose(did);

if (noise > 0)
    message = 'Outputting exodus file with noisy data'
    exo_put('noisy_data.exo', exodb);
end
     
     
