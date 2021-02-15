%%%%%
% This matlab file is used to check results against
% the unit test of the Axial Internal force.
% See UnitTestCode/IwanRegression.h for more details
%%%%%
clear all;
params.type='iwan';
%params.data=[chi(i) phi_max(i) R(i) S(i)];
%params.data=[-.6 7.1e-5 4.9e7 8.8e6];
params.data=[-0.4 1.71428571428571e-05 2894743368.62755 3333333.3333333];

alpha=[0 .25 0 .5]; 
dt=1e-5;
x0=0;
v0=0;
disp = 0;
force = 0;
Iwobj=Iwan(params.data);
format long;
forces = zeros(100,1); 
for i=1:100,
    x0 = x0 + 1e-8;
    [Iwobj,force]=Iwobj.updateForce(x0,v0);
    Iwobj=Iwobj.updateStates;
    %% update states in integrator and distributed damping
    forces(i)=force;
end
forces
%% head home.
datestr(now)
