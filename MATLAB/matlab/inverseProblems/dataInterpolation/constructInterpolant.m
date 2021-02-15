function [ Fu, Fv, Fw ] = constructInterpolant(filename)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

load(filename);

Fu = scatteredInterpolant(x0, y0, z0, nvar01);
Fv =scatteredInterpolant(x0, y0, z0, nvar02);
Fw = scatteredInterpolant(x0, y0, z0, nvar03);


end %function

