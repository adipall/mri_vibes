% this script reads in velocity data from a matlab file and compares a single-stage butterworth filter with a multistage butterworth filter

% first do the filter all in one shot
passFrequency = 200;
interp_ts = 5.0e-7;
time_coarsening = 1000;
butterCoeff = 2.0*interp_ts*passFrequency;
[bcoeff,acoeff] = butter(3,butterCoeff);
load v80a.mat;
vel = gvar66;
vel_filt = filter(bcoeff,acoeff,vel);
time_coarse = time(1:time_coarsening:length(time));
vel_filt_coarse = vel_filt(1:time_coarsening:length(time));
vel_coarse = vel(1:time_coarsening:length(time));


plot(time,vel)
hold
plot(time_coarse,vel_filt_coarse,'-r')
%hold


% now filter in 3 steps
% step 1, filter down to 25000Hz and coarsen by factor of 10
passFrequency = 20000;
interp_ts = 5.0e-7;
butterCoeff = 2.0*interp_ts*passFrequency;
[bcoeff,acoeff] = butter(3,butterCoeff);
vel_step1_filt = filter(bcoeff,acoeff,vel);
time_coarse_step1 = time(1:10:length(time));
vel_step1_filt_coarse = vel_step1_filt(1:10:length(time));
%
%% step 2, filter down to 2500Hz and coarsen by factor of 10
passFrequency = 2000;
interp_ts = 10*5.0e-7; % multiply by previous coarsening 
butterCoeff = 2*interp_ts*passFrequency;
[bcoeff,acoeff] = butter(3,butterCoeff);
vel_step2_filt = filter(bcoeff,acoeff,vel_step1_filt_coarse);
vel_step2_filt_coarse = vel_step2_filt(1:10:length(time_coarse_step1));
time_coarse_step2 = time_coarse_step1(1:10:length(time_coarse_step1));
%
%% step 2, filter down to 250Hz and coarsen by factor of 10
passFrequency = 200;
interp_ts = 10*10*5.0e-7; % multiply by previous coarsening 
butterCoeff = 2*interp_ts*passFrequency;
[bcoeff,acoeff] = butter(3,butterCoeff);
vel_step3_filt = filter(bcoeff,acoeff,vel_step2_filt_coarse);
vel_step3_filt_coarse = vel_step3_filt(1:10:length(time_coarse_step2));
time_coarse_step3 = time_coarse_step2(1:10:length(time_coarse_step2));
%
plot(time_coarse_step3,vel_step3_filt_coarse,'k')
xlabel('Time (s)')
ylabel('Surface velocity')
legend('unfiltered data','1-step filter','3-step filter')
axis([5.25 5.30 -80 80])

%%%%%%%
