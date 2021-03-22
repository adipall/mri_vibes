%% MRI Acquisition ; 2/25/2021
clear;
close all
clc

% read in voltage file
fname = "acq0063.csv"; 
T = readtable(fname);                   %, 'HeaderLines',3); 
t = T{:,1};
V1 = T{:,2};
V2 = T{:,3};

% read in FFT file
fname2 = "acq0062.csv"; 
T = readtable(fname2); 
t2 = T{:,1};
F_1 = T{:,2};
F_2 = T{:,3};

shift = t(1); % changes relative time to a more comprehensible start from 0.0s
if shift >= 0
    shift = -1*shift;
    t(:) = t(:) + shift;
else
    t(:) = t(:) + abs(shift);
end

tspan = t(end)-t(1);
ns = length(V1);
srate = ns/tspan;

figure(1);
subplot(2,1,1)
plot(t,V1) % FFT done by matlab
title('Channel 1')
% ylim([0,4e-3])
% xlim([0,1e4])
subplot(2,1,2)
plot(t, V2) % FFT done by waveforms
title('Channel 2')
% ylim([0,4e-3])
% xlim([0,1e4]

% Sampling frequency
Fs = srate;
% Perform fft
F = fft(V1);
% Length of fft
L = length(V1);

P2 = abs(F/L); % normalize by length
% Take first half of fft (amplitude)
P1 = P2(1:L/2+1);
% account for symmetry
P1(2:end-1) = 2*P1(2:end-1);
% 
% Frequencies associated with FFT 
fFreqs = Fs*(0:(L/2))/L;
response = max(P1);
% plot, but still need to use conversion factor to convert to disp/accel
figure(2);
subplot(2,1,1)
plot(fFreqs,P1) % FFT done by matlab
title('matlab')
ylim([0,4e-3])
xlim([0,1e4])
subplot(2,1,2)
plot(fFreqs,F_1) % FFT done by waveforms
title('waveforms')
ylim([0,4e-3])
xlim([0,1e4])
% 
% fileName = sprintf('trial_%d.mat',trialNumber);
% save(fileName,'fFreqs','ffData','response');


