function [b,a] = sdofwt(fn,sr,z,itype)
% sdofwt.m
% function [b,a]=sdofwt(fn,sr,z,itype)
%
% This routine finds the filter weights for a single-
% degree-of-freedom resonator using a ramp invarient filter.
% The weights are used in conjunction with the function
% sdoffltr and shspec to calculate the shock response
% spectrum of an acceleration time history using a ramp
% invarient simulation.
%
% Reference: 1. Smallwood D.O.,"An Improved Recursive Formula
%               for Calculating Shock Response Spectra," The
%               Shock and Vibration Bulletin 51(2), pp 211-217,
%               May 1981.
%            2. Smallwood D.O.,"The Shock Spectrum at Low Frequencies,"
%               The Shock and Vibration Bulletin 56(1), Appendix A,
%               pp 285-287, Aug 1986.
%
% Input: fn= the natural frequency of the sdof system in Hz
%        sr= the sample rate in samples/sec
%        z = the fraction of critical damping of the sdof system, 0<=z<1
%        itype = >0 sets up a base acceleration, absolute acceleration
%                   response system
%                <0 sets up a base acceleration, relative displacement
%                  (expressed in equivalent static acceleration units)
%                  system
%
% Output: b= three element row vector containing [B0 B1 B2] from ref. 2
%         a=three element row vector containing [1 A1P2 A2M1] from ref. 2

% David O. Smallwood, Sandia National Laboratories, Albuquerque NM
% June 6, 1988
% Copywrite 1994 Sandia National Laboratories

a = zeros(1,3);
b = zeros(1,3);
%pi=4*atan(1);
w=2*pi*fn/sr;
a(1) = 1;

if w<.0001
%if w < .001     % it may be necessary to change this constant
                % depending on the precision of the computer
% use these coefficients when the frequency is small for both
% models
%    small freq
     a(2) = 2*z*w + w*w*(1 - 2*z*z);
     a(3) = -2*z*w + 2*z*z*w*w;
     if itype > 0
%        abs accel model small freq
         b(1) = z*w + (w*w)*( (1/6) - 2*z*z/3 );
         b(2) = 2*w*w*(1-z*z)/3;
         b(3) = -z*w + w*w*( (1/6) - 4*z*z/3 );
     else
%        rel disp small freq
         b(1) = w*w/6;
         b(2) = 2*w*w/3;
         b(3) = w*w/6;
     end
else
%    large freq
     sq = sqrt(1 - z*z);
     e = exp(-z*w);
     wd = w*sq;
     sp = e * sin(wd);
     fact = (2*z*z -1)*sp/sq;
     c = e * cos(wd);
%    a1p2 & a2m1 are the same for both models
%    a1p2=a1+2=a(2)    a2m1=a2-1=a(3)
     a(2) = 2*(1-c);
     a(3) = -1 + e*e;
     if itype > 0
%         large freq abs accel
          spwd = sp/wd;
          b(1) = 1 - spwd;
          b(2) = 2*(spwd - e*cos(wd) );
          b(3) = e*e - spwd;
     else
%         large freq rel disp
          b(1) = ( 2*z*(c-1) + fact + w )/w;
          b(2) = ( -2*c*w + 2*z*(1-e*e) - 2*fact )/w;
          b(3) = ( e*e*(w+2*z) - 2*z*c + fact )/w;
     end
end
