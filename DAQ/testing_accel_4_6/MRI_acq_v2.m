%% MRI Acquisition ; 2/25/2021
clear;
close all
clc

% read in voltage file
fname = "test_full_0001.csv"; 
T = readtable(fname);
t = T{:,1};
V1 = T{:,2};
V2 = T{:,3};

shift = t(1); % changes relative time to a more comprehensible start from 0.0s
t(:) = t(:) + abs(shift);

tspan = t(end)-t(1);
end_loc = find(t == t(end));
ns = length(V1);
srate = ns/tspan;

figure(1) % if you wanted to see raw signal
plot(t,[V1,V2])
hold on
t_end = t(end_loc);
xlim([0,t_end])
legend('Channel 1', 'Channel 2')
xlabel('Time (s)')
ylabel('Amplitude (V)')

% Sampling frequency
Fs = srate;
% Perform fft
F = fft(V1);
% Length of fft
L = length(V1);

normed = abs(F/L); % normalize by length
% Take first half of fft (amplitude)
F_shifted = normed(1:L/2+1);
% account for symmetry
F_shifted(2:end-1) = 2*F_shifted(2:end-1);
% 
% Frequencies associated with FFT 
fFreqs = Fs*(0:(L/2))/L;
response = max(F_shifted);
% plot, but still need to use conversion factor to convert to disp/accel
figure(2)
plot(fFreqs,F_shifted)
xlim([0,500])
xlabel('Frequencies (Hz)')
ylabel('Amplitude (V)')


