%% MRI Acquisition ; 2/25/2021
clear;
close all
clc

%% read in best tube voltage file
fname = "acq0057.csv"; 
T = readtable(fname);
t = T{:,1};
V1 = T{:,2};
V2 = T{:,3};

shift = t(1); % changes relative time to a more comprehensible start from 0.0s
t(:) = t(:) + abs(shift);

tspan = t(end)-t(1);
ns = length(V1);
srate = ns/tspan

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

figure(1) % plot signal amplitude
subplot(1,3,1)
plot(fFreqs,F_shifted)
xlim([0,500])
ylim([0,0.01])
xlabel('Frequencies (Hz)')
ylabel('Amplitude (V)')
title('Specimen Tube')

%% read in best table voltage file
fname = "acq0009.csv"; 
T = readtable(fname);
t = T{:,1};
V1 = T{:,2};
V2 = T{:,3};

shift = t(1); % changes relative time to a more comprehensible start from 0.0s
t(:) = t(:) + abs(shift);

tspan = t(end)-t(1);
ns = length(V1);
srate = ns/tspan

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

figure(1) % plot signal amplitude
subplot(1,3,2)
plot(fFreqs,F_shifted)
xlim([0,500])
ylim([0,0.01])
xlabel('Frequencies (Hz)')
ylabel('Amplitude (V)')
title('Table')

%%
% read in best bore voltage file
fname = "acq0066.csv"; 
T = readtable(fname);
t = T{:,1};
V1 = T{:,2};
V2 = T{:,3};

shift = t(1); % changes relative time to a more comprehensible start from 0.0s
t(:) = t(:) + abs(shift);

tspan = t(end)-t(1);
ns = length(V1);
srate = ns/tspan

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

figure(1) % plot signal amplitude
subplot(1,3,3)
plot(fFreqs,F_shifted)
xlim([0,500])
ylim([0,0.01])
xlabel('Frequencies (Hz)')
ylabel('Amplitude (V)')
title('Bore')




