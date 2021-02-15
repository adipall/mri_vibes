% this script reads in velocity data and compares matlab's filtfilt with a manual reverse filter, but with 1000 times coarsening.

load v80a.mat
passFrequency = 200;
interp_ts = 5.0e-7;
butterCoeff = 2.0*interp_ts*passFrequency;
[bcoeff,acoeff] = butter(3,butterCoeff);
vel = gvar66;
vel_filt = filter(bcoeff,acoeff,vel);
vel_filt_filt = filtfilt(bcoeff,acoeff,vel);

%plot(time,vel)
%hold
plot(time,vel_filt)
hold
plot(time,vel_filt_filt,'g')

% now decimate and then reverse filter
time_coarse = time(1:1000:length(time));
vel_filt_coarse = vel_filt(1:1000:length(time));
coarse_butterCoeff = 1000*butterCoeff;
[bcoeff,acoeff] = butter(3,coarse_butterCoeff);
vel_filt_coarse_flip = flipud(vel_filt_coarse);
vel_filt_filt_coarse = filter(bcoeff,acoeff,vel_filt_coarse_flip);
plot(time_coarse,flipud(vel_filt_filt_coarse),'-r')

xlabel('Time (s)')
ylabel('surface velocity')
legend('filtered data','filtfilt','filtered, decimated by 1000, then reverse filtered')
axis([5.3 5.33 -60 30])

