%% MRI Acquisition ; 2/25/2021
clear;
close all
clc

%% read in voltage file
fname = "test_two_0017.csv"; 

T = readtable(fname);
t = T{:,1};
V1 = T{:,2}; % x-axis
V2 = T{:,3}; % z-axis

shift = t(1); % changes relative time to a more comprehensible start from 0.0s
t(:) = t(:) + abs(shift);

tspan = t(end)-t(1);
ns = length(V1);
srate = ns/tspan;
%% Conversion Factors
cf = 0.328125; % V/g
to_accel = 9.82/0.32; % (m/s^2)/V -> multiply V amp by this to get accel amp
%%
figure(1)
sgtitle('Raw Data and Spectral Analysis','FontSize',16,'FontWeight','bold')
subplot(1,2,1) % plot raw signal
plot(t,[V1,V2])
grid on
t_end = t(end);
xlim([0,t_end])
legend('X-Axis', 'Z-Axis','location','best')
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
fFreqs = (Fs*(0:(L/2))/L)';
response = max(F_shifted);

% throw out noise at 0 freq
fFreqs(1) = 0;
F_shifted(1) = 0;
[pks, locs] = findpeaks(F_shifted,'MinPeakProminence',0.9E-3);

figure(1)
subplot(1,2,2)
semilogx(fFreqs,F_shifted)
hold on
grid on
semilogx(fFreqs(locs),F_shifted(locs),'o')
for i = 1:size(locs,1)
    str_txt{i} = append(num2str(fFreqs(locs(i))),' Hz')
    shift_text(i) = fFreqs(locs(i))/4 % adjusts text position given frequency location (can't be uniform)
end
text(fFreqs(locs)+shift_text',F_shifted(locs),str_txt)

xlabel('Frequencies (Hz)')
ylabel('Amplitude (V)')
set(findall(gcf,'-property','FontSize'),'FontSize',14)
xlim([0,500]) % modify this, if you think there may be data past 500 Hz

% function peaks = plotter(x,y)
% 
% 
%
% end



