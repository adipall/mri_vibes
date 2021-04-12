%% MRI Acquisition ; 2/25/2021
clear;
close all
clc

% read in best tube voltage file
fname = "acq0057.csv"; 
T = readtable(fname);
t = T{:,1};
V1 = T{:,2};
V2 = T{:,3};

shift = t(1); % changes relative time to a more comprehensible start from 0.0s
t(:) = t(:) + abs(shift);
end_loc = find(t == t(end));
t_end = t(end_loc);

figure(1) % plot signal amplitude
subplot(1,3,1)
plot(t,[V1,V2])
hold on
xlim([0,t_end])
legend('Z-axis', 'X-axis','location','best','fontsize',10,'interpreter','latex')
xlabel('Time (s)','fontsize',14,'interpreter','latex')
ylabel('Amplitude (V)','fontsize',14,'interpreter','latex')
title('Specimen Tube','fontsize',16)

% read in best table voltage file
fname = "acq0009.csv"; 
T = readtable(fname);
t = T{:,1};
V1 = T{:,2};
V2 = T{:,3};

shift = t(1); % changes relative time to a more comprehensible start from 0.0s
t(:) = t(:) + abs(shift);
end_loc = find(t == t(end));
t_end = t(end_loc);

subplot(1,3,2)
plot(t,[V1,V2])
xlim([0,t_end])
legend('Z-axis', 'X-axis','location','best','fontsize',10,'interpreter','latex')
xlabel('Time (s)','fontsize',14,'interpreter','latex')
ylabel('Amplitude (V)','fontsize',14,'interpreter','latex')
title('Table','fontsize',16)

% read in best bore voltage file
fname = "acq0066.csv"; 
T = readtable(fname);
t = T{:,1};
V1 = T{:,2};
V2 = T{:,3};

shift = t(1); % changes relative time to a more comprehensible start from 0.0s
t(:) = t(:) + abs(shift);
end_loc = find(t == t(end));
t_end = t(end_loc);

subplot(1,3,3)
plot(t,[V1,V2])
xlim([0,t_end])
legend('Z-axis', 'X-axis','location','best','fontsize',10,'interpreter','latex')
xlabel('Time (s)','fontsize',14,'interpreter','latex')
ylabel('Amplitude (V)','fontsize',14,'interpreter','latex')
title('Bore','fontsize',16)

%%

tspan = t(end)-t(1);
ns = length(V1);
srate = ns/tspan;

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
ylim([0,0.01])
xlabel('Frequencies (Hz)')
ylabel('Amplitude (V)')




