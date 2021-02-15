function [mn,sd,sk,kt,xs,tt] = run_stats02(t,x,w,lap)
%RUN_STATS02
%
% Function to estimate the "running" statistics of x(t), a sample
% of a non-stationary process. Computed are time-varying estimates
% of the mean, standard deviation, skewness, and kurtosis of x(t).
% We also calculate the shifted and scaled version of x, defined by
% xs=(x-mn)/sd. This should be much faster than run_stats01.m.
%
%USAGE
%
% [mn,sd,sk,kt,xs,tt]=run_stats02(t,x,w,lap)
%
% t   = n x 1 time vector
% x   = n x 1 data vector
% w   = width (in time points) of averaging window
% lap = number of overlap points within each window

% 9/29/2010  (created)

% use Tim's code to produce array xw; each column in xw has w points,
% overlapping with lap points; zeros are added to last column if needed
xw = t_window(x(:)',w,lap);

% compute statistics column-wise
mn = mean(xw);
sd = std(xw);
sk = skewness(xw);
kt = kurtosis(xw);

% corresponding time vector
tt = mean(t_window(t(:)',w,lap));

% shifted/scaled version of x
xi = interp1(t,x,tt);   % interpolate original data onto tt
xs = (xi - mn)./sd;     % shift and scale xi

% outputs
mn=mn(:);sd=sd(:);
sk=sk(:);kt=kt(:);
xs=xs(:);tt=tt(:);
