clc
clear all

% script for generating different types of loading functions

total_time = 5e-5; % in seconds
nsteps = 1000;
frequency = 10000;

times = [0: total_time/nsteps : total_time];

y = sin(times * 2 * pi * frequency);

fid = fopen('./load.txt','w');
for i=1:length(times)
    fprintf(fid, 'data  %g %g\n', times(i), y(i));
end
plot (times, y);

fclose(fid);