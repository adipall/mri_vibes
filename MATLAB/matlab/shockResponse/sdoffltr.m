function y = sdoffltr(b,a,x,zi)
% function y = sdoffltr(b,a,x,zi)
% function filters a time history x with a single-degree-
% of-freedom oscillator using a ramp invarient filter where
% the filter weights [b,a] have been defined by the function
% sdofwt (see sdofwt.m)
%
% reference: sdofwt.m
%
% input: b= a three element vector of filter weights defined
%           by the function sdofwt
%        a= three element vector of filter weights defined
%           by the function sdofwt
%        x= the input vector to be filtered
%        zi= optional intial filter state (see filter), default is [0 0];
%
% output: y= the filtered time history
%
% note: the intial conditions for values of x&y less than index 1
%       are assumed to be zero

% revised 6/9/97 to add optional argument zi
% David Smallwood Sandia National Laboratories Albuquerque NM
% June 6, 1988
% Copywrite 1994 Sandia National Laboratories


if nargin<4,  zi=[0 0];, end
if length(zi)==0, zi=[0 0];, end
aa=[a(1) a(2)-2 a(3)+1];
y = filter(b,aa,x,zi);

