clear; close all; clc;

% Define number of elements and displacements
numelems = [5600 10761 32075 86088 103530 256600 359552];
avgydisp = [0.000216803 0.000817022 0.00188975 0.00287234 0.00398347 0.0069589 0.00220632];
maxydisp = [0.00021684 0.000817137 0.00189047 0.00287273 0.00398375 0.00696077 0.00220781];

% Initialize plot variables
n = length(numelems)-1;
avgydiff = zeros(1,n);
maxydiff = zeros(1,n);

for k = 1:n
    avgydiff(k) = abs(avgydisp(k+1)-avgydisp(k)) / avgydisp(k);
    maxydiff(k) = abs(maxydisp(k+1)-maxydisp(k)) / maxydisp(k);
end

figure;clf;
subplot(2,1,1)
plot(numelems(1:end-1),avgydiff, 'ko-')
xlabel('Total Number of Elements', 'interpreter', 'latex')
ylabel('Normalized Z-Displacement', 'interpreter','latex')
grid on
subplot(2,1,2)
plot(numelems(1:end-1),maxydiff, 'ko-')
xlabel('Total Number of Elements', 'interpreter', 'latex')
ylabel('Normalized Z-Displacement', 'interpreter','latex')
grid on