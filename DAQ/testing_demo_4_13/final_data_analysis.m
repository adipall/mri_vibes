%% Code to View Acquired Data & Frequency Spectrum - 4/13/2021
clear;
close all
clc

%% MAIN
% 30 Hz testers
% fname = "NR_30Hz_001";
% fname = "b5_r65_30Hz_001";
% fname = "b10_r80_30Hz_001";

% 45 Hz testers
% fname = "NR_45Hz_001";
% fname = "b5_r65_45Hz_001";
% fname = "b10_r80_45Hz_001";

% 60 Hz testers
% fname = "NR_60Hz_001";
% fname = "b5_r65_60Hz_001";
fname = "b10_r80_60Hz_001";



% prominence for Oscillator-type files - should be safe val, but can tweak
prom = 0.01;

run_funcs(fname,prom)

%% Functions
function run_funcs(fname,prom)
% Read data and do FFT calcs
read_data(fname)
% Run plotting code
plotter(prom,fname)
end

function read_data(fname)
fname_read = append(fname,".csv");
T = readtable(fname_read);
t = T{:,1};
V1 = T{:,2}; % x-axis
V2 = T{:,3}; % z-axis

shift = t(1); % changes relative time to a more comprehensible start from 0.0s
t(:) = t(:) + abs(shift);

tspan = t(end)-t(1);
ns = length(V2);
srate = ns/tspan;

%% FFT
Fs = srate; % Sampling frequency
% Perform fft
F = fft(V2);
% Length of fft
L = length(V2);

normed = abs(F/L); % normalize by length
% Take first half of fft (amplitude)
F_shifted = normed(1:L/2+1);
% account for symmetry
F_shifted(2:end-1) = 2*F_shifted(2:end-1);

% Frequencies associated with FFT 
fFreqs = (Fs*(0:(L/2))/L)';
response = max(F_shifted);

% throw out noise at 0 freq
fFreqs(1) = 0;
F_shifted(1) = 0;

save('for_plotting.mat')
end

function plotter(prom,fname) 
    load('for_plotting.mat')
    figure(1)
    sgtitle(fname,'FontSize',16,'FontWeight','bold','Interpreter', 'none')
    subplot(1,2,1) % plot raw signal
    plot(t,[V1,V2])
    grid on
    t_end = t(end);
    xlim([0,t_end])
    legend('X-Axis', 'Z-Axis','location','best')
    xlabel('Time (s)')
    ylabel('Amplitude (V)')
    
    [pks, locs] = findpeaks(F_shifted,'MinPeakProminence',prom); % find peaks using specified prominence
    
    figure(1)
    subplot(1,2,2)
    semilogx(fFreqs,F_shifted) % plot spectrum
    hold on
    grid on
    semilogx(fFreqs(locs),F_shifted(locs),'o') % plot peaks
    
    for i = 1:size(locs,1)
        str_txt{i} = append(num2str(fFreqs(locs(i))),' Hz; ', num2str(F_shifted(locs(i))),' V')
        shift_text(i) = fFreqs(locs(i))/4 % adjusts text position given frequency location (can't be uniform)
    end
    text(fFreqs(locs)+shift_text',F_shifted(locs),str_txt)

    xlabel('Frequencies (Hz)')
    ylabel('Z-Axis Amplitude (V)')
    set(findall(gcf,'-property','FontSize'),'FontSize',14)
    xlim([0,500]) % can modify this if want better x-resolution
    set(gcf,'position',[0 100 1000 800])
    savename = append(fname,".png");
    saveas(gcf,savename)
end