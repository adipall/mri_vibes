function [f,k,c,y0]=iwan_function_preferred(x,v,y0,params)
%
%  INPUT PARAMETERS:

% params = [ chi, beta, K_0, F_S];
% This is the preferred parameters set.
% All terms are dimensionless or of integral dimension.

% For the purpose of this calculation parameter set
% [ chi, phi_max, R, S] is more useful and will be calculted


%   y0 are state variable.
%      = slider positions (nsliders+1) (+1 for the macroslip)
%   Before  first call to this routine set            y0=zeros(30,1);

%  x is relative displacement imposed on Iwan element (a number)
%  v is relative velocity imposed on element.  v is not used here

%#ok<*INUSL>

% OUTPUT PARAMETERS
%  f = force in element
%  k = current value of tangent stiffness
%  c = current viscous damping  == 0
% y0, current value of state variables for use at next time step.

chi = params(1);
beta= params(2);
K_T = params(3);
F_S = params(4);

% evaluate terms of other parameter set
phi_max = F_S*(1+beta)/(K_T*(beta+(chi+1)/(chi+2)));
R = F_S*(chi+1)/((phi_max^(chi+2))*(beta+(chi+1)/(chi+2)));
S = (F_S/phi_max)*beta/(beta+(chi+1)/(chi+2));


% Written by Todd Simmermacher 10/2004
% modified by DJS 12 November 2005 for generic function call
% comments modified again by DJS 2 July 2007
% Changed to take preferred parameter set by DJS 18 June 2009
% Updated by MRWB February 2011

c=0;

n_sliders=length(y0);
alpha=1.2;

%ju=linspace(0,phi_max,n_sliders); % for a linear spacing of integration
%points

% the following sets up the non uniform integration points
dphi=phi_max*(alpha-1)/(alpha^(n_sliders-1) -1);
ju=cumsum([0 dphi*1.2.^(0:n_sliders-2)]);
phi=ju(1:end-1)+diff(ju)/2;
phi=[phi(:);phi_max];

% define some commonly used quantities
dif=x-y0;
phi1=ju(:).^(chi+1);
phi2=ju(:).^(chi+2);
dif=dif(:);

% figure out what sliders are sliding
Il=find((abs(dif(1:end-1))-phi(1:end-1))<0); %fabs(u-x)<phi
Ig=find((abs(dif(1:end-1))-phi(1:end-1))>=0); %fabs(u-x)>=phi

%if nargout==3,  % calculate tangent stiffness
    k=0;
    if ~isempty(Il),
        f1=R*sum((phi1(Il+1) - phi1(Il)))/(chi+1);
        if abs(dif(end))<phi_max,
            k=f1+S;
        else
            k=f1;
        end
    end
%end

% some more useful quantities
Iz= ~dif;
dif(Iz)=1e-30;
dir=dif./abs(dif);

clear f

% calculate the force in the element
f1=R*sum(dif(Il).*(phi1(Il+1)-phi1(Il)))/(chi+1); %not sliding contribution
f2=R*sum(dir(Ig).*(phi2(Ig+1)-phi2(Ig)))/(chi+2); %sliding contribution
f=f1+f2;

% add in the macroslip singularity
if abs(dif(end))<phi(end),
    f=f+S*dif(end);
else
    f=f+dir(end)*S*phi_max;
end

% update states
y0(Ig)=x-dir(Ig).*phi(Ig);
% the next line contains the correction to bug discoverd by Wil
%                    Holzman and Todd Simmermacher
if abs(dif(end)) > phi(end),
    y0(end)=x-dir(end)*phi(end);
end



