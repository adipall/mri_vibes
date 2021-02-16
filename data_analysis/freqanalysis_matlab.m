clear;
close all
clc

load('trial_2.mat'); % manually change to appropriate name of file
%load('controlDataex_3.mat')
% constants and independent variables
len = length(response);         % for indexing and freq vector
maxfreq = 3200;
minfreq = 1200;
inc = 50;
freqs = minfreq:inc:maxfreq;          % freq values
srate = 6400;                   % samplingrate
tfull = (0:(1/srate):10)';      % orignal time of signal
htloc = length(tfull)/2;        % index of half time
t = tfull(1:htloc);             % half time of signal
cf = 1.02e-4;                   % conversion factor V/(m/s^2)

% output set-up
a_amp = zeros(len,1);
u_amp = zeros(len,1);

for f= 1:len
    w = freqs(f)*2*pi; 
    V = real(response(f)*exp(-1i*w*t));
    a_amp(f) = response(f)/cf; % response(f) is the voltage amplitude
    u_amp(f) = a_amp(f)/(w^2);
    
end
freqs = freqs'; % transposed freqs to keep rows(55)consistent with all other matrices 

% manually change analysisresultsn.mat
response = response';
data = [freqs response u_amp a_amp];
save('mockdata.mat','freqs','response','u_amp','a_amp'); 
csvwrite('mockdata.csv',data)
