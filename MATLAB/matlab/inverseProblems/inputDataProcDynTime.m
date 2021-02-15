clc
clear all

inputString = './data';
truthtablefile = './ttable.txt';
datafile = './data.txt';
% Measured dofs
X = 1;
Y = 1;
Z = 1;
locations = 9;
noise = 0.0;


inputfile = [inputString,'0.csv'];
importfile(inputfile);

%data = load(inputfile);

nodes(1:locations) = 0;
nodes(1) = data(1,4);
dt = data(2,6) - data(1,6);
steps = length(data(:,1));

%Construct the data file
did = fopen(datafile,'w');
  
fprintf(did, '%d  %d %d\n', 3*locations, steps, dt);
   
temp(1:3*locations, steps)=0;

for p=1:locations
    inputfile = [inputString,int2str(p-1),'.csv'];
    importfile(inputfile);
    %data = load(inputfile);
           
    dispx = data(:,1) .* (1 + noise*randn(steps,1));   
    dispy = data(:,2) .* (1 + noise*randn(steps,1));
    dispz = data(:,3) .* (1 + noise*randn(steps,1));
    nodes(p) = data(1,4);
        
    for s=1:steps  
        temp(p,s) = dispx(s);
        temp(p+locations,s) = dispy(s);
        temp(p+2*locations, s) = dispz(s);
    end   
end

for j=1:3*locations
    for p=1:steps
        fprintf(did, '%g ', temp(j,p));
    end
    fprintf(did,'\n');
end

%Construct the truth table
tid = fopen(truthtablefile,'w');
fprintf(tid, '%d\n', locations);

for i=1:locations
    fprintf(tid,'%g ', nodes(i));
    fprintf(tid, '%d %d %d\n', X, Y, Z);
end
fclose(tid);


fclose(did);
     
     
