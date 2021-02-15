clc
clear all

format long

%****************************INPUT PARAMETERS******************
inputfile ='./measured_data.e';
nodeset = 1;

noise = 0.0;
%Outputfile names
truthtablefile = './ttable.txt';
datafile = './dataReal.txt';
datafile_im = './dataImag.txt';
%***************************************************************
exodb = exo_get(inputfile);
nodes = exodb.Nodesets(nodeset).Nodes;
locations = length(nodes);
steps = length(exodb.Time);

%Construct the data file
did = fopen(datafile,'w');
did2 = fopen(datafile_im,'w');
  
fprintf(did, '%d  %d\n', locations, steps);
fprintf(did2, '%d  %d\n', locations, steps);

node_num_map = exodb.Nodes.NodeNumMap;
num_nodes_in_mesh = length(exodb.Nodes.Coordinates);
if length(node_num_map)==0
   for i=1:num_nodes_in_mesh
      node_num_map(i) = i;
   end
end



    for p=1:locations
        for s=1:steps
            dispr = exodb.NodalVars(1).Data(nodes(p), s)*(1+noise*rand());
            dispi = exodb.NodalVars(2).Data(nodes(p), s)*(1+noise*rand());      
            fprintf(did, '%17.12e     ', dispr);
            fprintf(did2, '%17.12e     ', dispi);
        end
        fprintf(did,'\n');
        fprintf(did2,'\n');
    end



%Construct the truth table
tid = fopen(truthtablefile,'w');
fprintf(tid, '%d\n', locations);

for i=1:locations
    mapped_node = node_num_map(nodes(i))
    fprintf(tid,'%d\n', mapped_node);
end
fclose(tid);

fclose(did);
     
     
