function inputDataProcDynFreq_Exo(nodesetId)
%addpath(genpath('/data1/opt/sierra/matlab/'));
if (nargin < 1)
   nodesetId=1; 
end

clc
%clear all

format long

%****************************INPUT PARAMETERS******************
inputfile ='simplemodel-dfrf.e';
% Measured dofs
X = 1;
Y = 1;
Z = 1;
noise = 0.0;
%Outputfile names
truthtablefile = './ttable.txt';
datafile = './dataReal.txt';
datafile_im = './dataImag.txt';
%***************************************************************
exodb = exo_get(inputfile);

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

%Construct the data file
did = fopen(datafile,'w');
did2 = fopen(datafile_im,'w');
  
fprintf(did, '%d  %d\n', 3*locations, steps);
fprintf(did2, '%d  %d\n', 3*locations, steps);

node_num_map = exodb.Nodes.NodeNumMap;
num_nodes_in_mesh = length(exodb.Nodes.Coordinates);
if length(node_num_map)==0
   for i=1:num_nodes_in_mesh
      node_num_map(i) = i;
   end
end

for d=1:3
    for p=1:locations
        for s=1:steps
            dispr = exodb.NodalVars(d).Data(nodes(p), s)*(1+noise*randn());
            dispi = exodb.NodalVars(d+3).Data(nodes(p), s)*(1+noise*randn());      
            fprintf(did, '%17.12e     ', dispr);
            fprintf(did2, '%17.12e     ', dispi);
        end
        fprintf(did,'\n');
        fprintf(did2,'\n');
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
     
     
