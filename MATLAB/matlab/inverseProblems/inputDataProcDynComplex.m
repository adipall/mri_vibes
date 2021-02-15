clc
clear all

inputString = './data3.';
truthtablefile = './ttable.txt';
datafile = './data.txt';
datafile_im = './data_im.txt';
tsteps = 1;
noise = 0.0;


inputfile = [inputString,'0.csv'];
data = load(inputfile);

nodes = data (:,7);
no_nodes = length(nodes);

%Construct the truth table
tid = fopen(truthtablefile,'w');
fprintf(tid, '%d\n', no_nodes);

for i=1:no_nodes
    fprintf(tid,'%g ', nodes(i));
    fprintf(tid, '%d %d %d\n', 1, 1, 1);
end
fclose(tid);




%Construct the data file
did = fopen(datafile,'w');
did2 = fopen(datafile_im,'w');
  
fprintf(did, '%d  %d\n', 3*no_nodes, tsteps);
fprintf(did2, '%d  %d\n', 3*no_nodes, tsteps);
   
temp(1:3*no_nodes, tsteps)=0;
tempim(1:3*no_nodes, tsteps)=0;

for s=1:tsteps
    inputfile = [inputString,int2str(s-1),'.csv'];
    data = load(inputfile);
    
    dispx = data(:,1) .* (1 + noise*randn(no_nodes,1));   
    dispy = data(:,3) .* (1 + noise*randn(no_nodes,1));
    dispz = data(:,5) .* (1 + noise*randn(no_nodes,1));
    dispImx = data(:,2) .* (1 + noise*randn(no_nodes,1));
    dispImy = data(:,4) .* (1 + noise*randn(no_nodes,1));
    dispImz = data(:,6) .* (1 + noise*randn(no_nodes,1));
    

    for p=1:no_nodes
        temp(p,s) = dispx(p);
        temp(p+no_nodes,s) = dispy(p);
        temp(p+2*no_nodes, s) = dispz(p);
        tempim(p,s) = dispImx(p);
        tempim(p+no_nodes,s) = dispImy(p);
        tempim(p+2*no_nodes, s) = dispImz(p);
    end   
end

for j=1:3*no_nodes
    for p=1:tsteps
        fprintf(did, '%g ', temp(j,p));
        fprintf(did2, '%g ', tempim(j,p));
    end
    fprintf(did,'\n');
      fprintf(did2,'\n');  
end

fclose(did);
     
     
