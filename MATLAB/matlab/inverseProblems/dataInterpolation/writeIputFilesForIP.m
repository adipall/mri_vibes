function [] = writeIputFilesForIP(nodesWithData,node_num_map, dispX, dispY, dispZ)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Outputfile names

truthtablefile = './ttable.txt';
datafile = './data.txt';
datafile_im = './data_im.txt';
X=1; %Set this to zero to turn off the dof
Y=1;
Z=1;
%********************************
locations = length(nodesWithData);
steps=1;

disps=[dispX dispY dispZ];

%Construct the data file
did = fopen(datafile,'w');
did2 = fopen(datafile_im,'w');
  
fprintf(did, '%d  %d\n', 3*locations, steps);
fprintf(did2, '%d  %d\n', 3*locations, steps);

for d=1:3
    for p=1:locations
        dispr = disps(p,d);
        dispi = 0;      
        fprintf(did, '%17.12e     ', dispr);
        fprintf(did2, '%17.12e     ', dispi);
        fprintf(did,'\n');
        fprintf(did2,'\n');
    end
end


%Construct the truth table
tid = fopen(truthtablefile,'w');
fprintf(tid, '%d\n', locations);

for i=1:locations
    mappedNode=node_num_map(nodesWithData(i));
    fprintf(tid,'%d  ', mappedNode);
    fprintf(tid, '%d %d %d\n', X, Y, Z);
end
fclose(tid);
fclose(did);
fclose(did2);
end

