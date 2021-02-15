function [x,e]=modaltran(dt,nsteps,modalforce,timefunc,freq,damp,nvars,evars)
% [x,e]=modaltran(dt,nsteps,modalforce,timefunc,freq,damp, nvars, evars) 
%
% computes the time response of a structure using the modal
% decomposition. This is designed to interface nicely with
% Salinas.
%
%
% From the Salinas history file we get the following:
%   nvar01, nvar02, etc. These are the nodal values at the points of interest.
%   evar01, evar02, etc. These are the element results for the modes
%   EigenFrequency - frequency (in Hz) of solution. Usually gvar02
%
%   nvars=[nvar01' nvar02' ...]';
%   evars=[evar01' evar02' ...]';
%
% note that some element variables (such as vonmises stress)
% should not be linearly combined. This routine blindly does
% such a combination.
%
% dt is the time step
% nsteps is the number of steps to integrate
% modalforce(nmodes,nforces) is a matrix of modal forces
% timefunc(nforces,1) is an array of time functions
% freq. modal frequency in Hz (nmodes,1)
% damp(nmodes,1) -  modal damping coefficient (like 0.02 for 2% damping)
% nvars, evars described above
%
% modalforce is available as a separate m file from salinas. It is 
% the modal force vector for a unit time function.


% $Revision$
% $Date$

Lambda = (2*pi*freq).^2;
delh=dt/2.;
delsq=delh^2;
nmodes=size(modalforce,1);
cdamp=damp.*freq*2*pi*2;

qforce=zeros(nmodes,1);
q_n=zeros(nmodes,1);
qd_n=zeros(nmodes,1);
qdd_n=zeros(nmodes,1);
xdim=size(nvars,1);
x=zeros(xdim,nsteps);

for it=1:nsteps
  t=it*dt;
  tf=eval(sprintf('%s(%g)',timefunc,t));
  qforce=modalforce*tf;

  qf2 = qforce - cdamp.*(qd_n+delh*qdd_n) - ...
      Lambda.*(q_n+(2*delh*qd_n)+delsq*qdd_n);
  
  qddplus = qf2./( 1 + delh*cdamp+delsq*Lambda);
  
  q_n = q_n + 2*delh*qd_n + delsq*(qdd_n+qddplus);
  qd_n = qd_n + delh*(qdd_n+qddplus);
  qdd_n = qddplus;
  x(:,it)=nvars*q_n;
  e(:,it)=evars*q_n;
end
