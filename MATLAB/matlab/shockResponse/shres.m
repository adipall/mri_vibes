function x = shres(a,t,fn,ZETA)
% shres.m
% function x = shres(a,t,fn,ZETA)
%
% function to calculate the residual shock response spectrum
%   The residual response is the peak response of the sdof system
%   as it decays after the input has ended, i.e., the input is zero.
%   The first two peaks in the decaying response will be the largest.
%   One peak will be negative and one positive.  We don't know before
%   the calculation if the first and largest peak in amplitude will be
%   positive or negative.
% INPUT: a = a two element vector giving the response at two times
%            after the input is zero
%        t = a two element vector giving the times where a was evaluated
%        fn = the natural frequency, inverse units of time, Hz if time is
%             sec.
%        ZETA =  the fraction of critical damping for the sdof filter
% OUTPUT:
%        x = a two element vector giving the amplitude of the response at
%            the first two peaks
% Method:
% The response has the form
% a(t) = exp(-ZETA wn t)[z(1)sin(wd t) + z(2)cos(wd t)].
% If I know a(t) at two values of time, t, I can calculate the constants
% z(1) and z(2).  The general form of the response can then be solved for
% the maximum by:  Finding the time of the first maximum by setting the
% derivative to zero.  Then substituting the time of the maximum response
% back into the general equation to find the maximum response.
% The second peak will occur pi radians later.

% D. Smallwood, Sandia National Laboratories, August 4, 1989
% Copywrite 1994 Sandia National Laboratories
% modified 10/24/96 to fix minor bug in the calculation of tmax;
%     tmax must be between 0 and pi/wd.
% modified 5/13/97 to intialize x.
% modifed 7/8/98 to return zero if input is zero
% modified 8/24/98 to correct some of the error checks, and add warning

   x = [0; 0];
if ZETA>1
   error(['shres: ZETA>1: ZETA = ' num2str(ZETA)])
   x = [0; 0];
   return
end
if ZETA<0
   error(['shres: ZETA<0: ZETA = ' num2str(ZETA)])
   x = [0; 0];
   return
end
if fn<0
   error(['shres: fn<0: fn = ' num2str(fn)])
   x = [0 ; 0];
   return
end

ac = a(:);                    % make sure ac is a column vector
if max(abs(ac))<=eps
   x = [0; 0];
   return
end

if abs(t(2)-t(1))>=1/fn
   x = a;
   warning(['time difference larger than a period of fn: fn=' num2str(fn) '  t(2)-t(1)=' num2str(t(2)-t(1))])
   return
end
pi2 = 2 * pi;
fd = fn * sqrt(1-ZETA*ZETA);  % find the damped natural frequency
wd = pi2 * fd;                % damped freq in radians/sec
wn = pi2 * fn;                % natural freq in radians/sec
wdt1 = wd*t(1);
wdt2 = wd*t(2);
e1 = exp(-wn*ZETA*t(1));
e2 = exp(-wn*ZETA*t(2));
c = zeros(2,2);
sd1 = sin(wdt1);
sd2 = sin(wdt2);
cd1 = cos(wdt1);
cd2 = cos(wdt2);
c(1,1) = e1*sd1;
c(1,2) = e1*cd1;
c(2,1) = e2*sd2;
c(2,2) = e2*cd2;
% make sure x and t are a column vectors
xx = x(:);
tt = t(:);
z = c\ac;
piwd = pi/wd;

% the response has the form
% a(t) = exp(-ZETA wn t)[z(1)sin(wd t) + z(2)cos(wd t)]
% the first maximum response will occur at the time, tmax
% tmax = (1/wd) * arctan((wd z(1) - ZETA wn z(2))/(wd z(2) - ZETA wn z(1)))

tmax = (1/wd) * atan( (wd*z(1) -ZETA*wn*z(2))/(wd*z(2) -ZETA*wn*z(1)) );
if tmax<0
   tmax = tmax + piwd;
end
if tmax>piwd
   tmax = tmax - piwd;
end
%  find the maximum

x(1) = exp(-ZETA*wn*tmax)*( z(1)*sin(wd*tmax) + z(2)*cos(wd*tmax) );

% the next peak will occur pi radians later in time

tmin = tmax + piwd;
x(2) = exp(-ZETA*wn*tmin)*( z(1)*sin(wd*tmin) + z(2)*cos(wd*tmin) );

