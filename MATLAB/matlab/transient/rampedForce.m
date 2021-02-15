%RAMPEDFORCE function f=
% m u'' + c u' + k u = vkt
% a time integrator has a function to evaluate the
% force at a time.  v is a parameter.
function f= rampedForce(k,v,tf)
f=v*k*tf;
