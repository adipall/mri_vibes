clear; close all; clc;

% Define number of elements and displacements
numelems = [16606 17050 19763 23718 33152 51961 130856];
avgydisp = [1 2 25 3 35 4 5];
maxydisp = [1 2 25 3 35 4 5];

% Initialize plot variables
n = length(numelems)-1;
avgydiff = zeros(1,n);
maxydiff = zeros(1,n);

for k = 1:n
    avgydiff(k) = abs(avgydisp(k+1)-avgydisp(k));
    maxydiff(k) = abs(maxydisp(k+1)-maxydisp(k));
end

figure;clf;
subplot(2,1,1)
plot(numelems(1:end-1),avgydiff)
subplot(2,1,2)
plot(numelems(1:end-1),maxydiff)
