function [s,fn]=shspec(x,fn,z,sr,itype,bwt,awt)
% shspec.m
% function [s,fn] = shspec(x,fn,z,sr,itype,bwt,awt)
%
% Function to calculate the shock response spectrum using a ramp 
% invarient digital filter simulation of the single-degree-of-freedom
% -system.
%
% INPUT:   x = a vector containing the input time history
%         fn = a vector containing the frequencies at which the spectrum
%              will be calculated.
%              Default: fn=logspace(log10(sr/1e4),log10(sr/4),50);
%          z = the fraction of critical damping for the spectrum
%              calculation. Default: z=.05
%         sr = the sample rate of the time history, x. Default: sr=1
%      itype = the type of spectrum desired
%            > 0 (pos)-- base acceleration-absolute acceleration model
%            < 0 (neg)-- base acceleration-relative displacement model
%                   (expressed in equivalent static acceleration units)
%            If abs(itype) is:
%            1--positive primary,  2--negative primary,  3--absolute maximum primary
%            4--positive residual, 5--negative residual, 6--absolute maximum residual
%            7--largest of 1&4, maximum positive, 8--largest of 2&5, maximum negative
%            9 -- maximax, the largest absolute value of 1-8
%           10 -- returns a matrix s(9,length(fn)) with all the types 1-9.
%            Default: itype=9
%      bwt,awt = optional weights length(fn)x3. If not given, computed by sdofwt.m 
% All input arguments except the first are optional.
% If an input is empty or missing the default value is used.
%
%  OUTPUT: s = a matrix containing the shock response spectrum at the
%              frequencies, fn.  If abs(itype)<10, s will have 1 row
%              If abs(itype)==10 each row of s will contain 1 spectrum
%              type for all freq's; each column will contain all types for one freq.
%          fn = the natural frequencies at which the shock spectrum is calculated.

%
%  other routines used: sdofwt -- designs the filters
%                       sdoffltr -- filters a time history with a sdof
%                                   filter
%                       shres -- finds the residual response
%
% Reference: 1. Smallwood D.O.,"An Improved Recursive Formula
%               for Calculating Shock Response Spectra," The
%               Shock and Vibration Bulletin 51(2), pp 211-217,
%               May 1981.
%            2. Smallwood D.O.,"The Shock Spectrum at Low Frequencies,"
%               The Shock and Vibration Bulletin 56(1), Appendix A,
%               pp 285-287, Aug 1986.
%
% D. O. Smallwood, Sandia National Laboratories, 19 June 1989.
% Modified 8-4-89 to add itype=10
% Added optional input weights bwt,awt 3/4/94 dos
% Added error checking, defaults, and 2nd output variable, 6/8/98, dos
% Copywrite 1994-1998 Sandia National Laboratories

% set defaults
if nargin<5, itype=9;, end
if isempty(itype), itype=9;, end
if nargin<4, sr=1;, end
if isempty(sr), sr=1;, end
if nargin<3, z=.05;, end
if isempty(z), z=.05;, end
if nargin<2, fn=logspace(log10(sr/1e4),log10(sr/4),50);, end
if isempty(fn),fn=logspace(log10(sr/1e4),log10(sr/4),50);, end
   
nf=length(fn);
nx=length(x);
if nf==0
   disp('Error in shspec: length of vector of natural frequencies, fn==0')
   return
end
if nx==0
   disp('Error in shspec: length of input vector x is zero')
   return
end
xmax = max(abs(x));
jtype=abs(itype);
   if xmax<=eps        % special case when input is zero
      if jtype < 10, s(1,:) = zeros(1,nf);, end
      if jtype == 10, s = zeros(9,nf);, end
      return
   end

xx = x(:)';                   % make sure time history is a row vector
errorshspec = 'error in shspec';
if jtype>10, error([errorshspec 'itype= ' num2str(itype)]), end
if itype==0, error([errorspec ' itype is zero']), end
if sr<=0, error([errorshspec ' sr<=0, sr= ' num2str(sr)]), end
if min(fn)<=0, error, fn, end
if z<0, error([errorshpec 'damping is negative z= ' num2str(z)]), end 
if z==0, warning('shspec: damping is zero'), end
if jtype<10, ns=1; end
if jtype==10, ns=9; end
s=zeros(ns,nf);

% loop for each natural frequency
for i=1:nf
    if nargin<7
       [b,a] = sdofwt(fn(i),sr,z,itype);
    elseif isempty(bwt)|isempty(awt)
       [b,a] = sdofwt(fn(i),sr,z,itype); 
    else
       b = bwt(i,:);, a = awt(i,:);
    end
%   append enough zeros to get 1/100 of the residual response
    nzeros = max(2, ceil(.01*sr/fn(i)) );
    xe = [xx zeros(1,nzeros)];

%   filter the appended input to find the filter output at one natural frequency
    y =  sdoffltr(b,a,xe);
%disp('shspec: line 102'), keyboard

%   find primary response

    primax=max(y(1:nx));
    primax=max(0,primax);
    primin=min(y(1:nx));
    primin=abs(min(0,primin));
    priabs=max(primax,primin);

%   find residual response

    t(1)=0;
    t(2)=(nzeros-1)/sr;
    zy=zeros(1,2);
    zy(1)=y(nx+1);
    zy(2)=y(nx+nzeros);
    res=zeros(1,2);
    %if abs(t(2)-t(1))>=1/fn(i), error('shres: time difference too large for fn'), end
    res = shres(zy,t,fn(i),z);
%disp('shspec: line 120'), keyboard
    respos=max(res);
    respos=max(0,respos);
    resneg=min(res);
    resneg=abs(min(0,resneg));
    resabs=max(respos,resneg);
    ss(1) = primax;
    ss(2) = primin;
    ss(3) = priabs;
    ss(4) = respos;
    ss(5) = resneg;
    ss(6) = resabs;
    ss(7) = max(respos,primax);
    ss(8) = max(resneg,primin);
    ss(9) = max(priabs,resabs);
    if jtype < 10, s(1,i)=ss(jtype); end
    if jtype == 10, s(:,i)=ss'; end
end

