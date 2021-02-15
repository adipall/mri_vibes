function [force, RIPP] = RIPPjoint(d,v,RIPP)
% This function calculates the force acting across a bolted joint. The
% formulation for this model is based on a Reduced Iwan model that includes
% Pinning (Reduced Iwan Plus Pinning: RIPP). This function takes as input:
%    d : instantaneous displacement across the joint
%    v : instantaneous velocity across the joint
% RIPP : a structure that includes
%          FS : Force necessary to induce macroslip
%          KT : Tangent stiffness in the microslip regime
%         chi : Dissipation behavior of the joint ( -1 < chi <= 0)
%        beta : Ratio of sliders that slip early compared to late
%       theta : Ratio of dynamic to static friction (0 < theta <= 1)
%          Kp : Pinning stiffness
%          dp : Pinning engagement distance
%       model : 1 - Segalman/Mignolet's model,
%               2 - a uniform distribution,
%               3 - Iwan-Stribeck 8-parameter model
%                     This requires Fv (viscous damping due to lubrication)
%                                   vs (scaling parameter for transition)
%                                   ds (exponent for transition)
%      PhiMax : Maximum displacement for macroslip. This is not an
%               independent variable, and is initialized via
% RIPP.PhiMax = RIPP.FS*(1+RIPP.beta)/RIPP.KT/(RIPP.beta+(RIPP.chi+1) ...
%               /(RIPP.chi+2)) ;
% % Other initialization variables can be setup via the script
% RIPP.F0 = 0 ;
% RIPP.d0 = 0 ;
% RIPP.direction = 0 ;
% RIPP.Flast = 0 ;
% RIPP.dlast = 0 ;
%
% % A sample set of parameters might be
% RIPP.Kp = 2e7 ;
% RIPP.dp = 2e-3 ;
% RIPP.FS = 4e3 ;
% RIPP.KT = 1.5e7 ;
% RIPP.chi = -0.5 ;
% RIPP.beta = 0.005 ;
% RIPP.theta = 0.8 ;
% RIPP.PhiMax = RIPP.FS*(1+RIPP.beta)/RIPP.KT/(RIPP.beta+(RIPP.chi+1)/(RIPP.chi+2)) ;
% % And for the Iwan-Stribeck model
% RIPP.Fv = 1e-3 ;
% RIPP.ds = 2 ;
% RIPP.vs = 1e-4 ;

% Determine if there is a direction change
if sign(v) ~= RIPP.direction
  RIPP.F0 = RIPP.Flast ;
  RIPP.d0 = RIPP.dlast ;
  RIPP.direction = sign(v) ;
end
RIPP.dlast = d ;

% Define the relative coordinate system and scale appropriately
u = min(sign(v)*(d - RIPP.d0)*RIPP.FS/(RIPP.FS-sign(v)*RIPP.F0),RIPP.PhiMax) ;

% Calculate pinning forces
Fpin = 0 ;
if d > RIPP.dp
  Fpin = RIPP.Kp*(d-RIPP.dp) ;
elseif d < -RIPP.dp
  Fpin = RIPP.Kp*(d+RIPP.dp) ;
end

% Calculate Iwan forces
if RIPP.model == 1 % Segalman 2005 (theta = 1), Mignolet 2015 (theta != 1)
  FIwan = RIPP.FS*(RIPP.chi+1) ...
         /(RIPP.PhiMax^(RIPP.chi+2)*(RIPP.beta+(RIPP.chi+1)/(RIPP.chi+2))) ...
         *((RIPP.theta/(RIPP.chi+2)-1/(RIPP.chi+1))*u^(RIPP.chi+2) ...
          +RIPP.PhiMax^(RIPP.chi+1)/(RIPP.chi+1)*u) ...
         +RIPP.FS/RIPP.PhiMax*RIPP.beta/(RIPP.beta+(RIPP.chi+1)/(RIPP.chi+2)) ...
         *min(u,RIPP.PhiMax) ;
elseif RIPP.model == 2 % Uniform Iwan distribution
  FIwan = 4/3*RIPP.FS/RIPP.PhiMax*(u-1/4/RIPP.PhiMax*u^2) ;
elseif RIPP.model == 3 % Iwan-Stribeck model
  FIwan = RIPP.FS*(RIPP.chi+1) ...
          /(RIPP.PhiMax^(RIPP.chi+2)*(RIPP.beta+(RIPP.chi+1)/(RIPP.chi+2))) ...
          *((RIPP.theta/(RIPP.chi+2)-1/(RIPP.chi+1))*u^(RIPP.chi+2) ...
          +RIPP.PhiMax^(RIPP.chi+1)/(RIPP.chi+1)*u) ...
        +RIPP.FS/RIPP.PhiMax*RIPP.beta/(RIPP.beta+(RIPP.chi+1)/(RIPP.chi+2)) ...
          *min(u,RIPP.PhiMax) ...
        +1/RIPP.PhiMax^(RIPP.chi+1)/(RIPP.beta+(RIPP.chi+1)/(RIPP.chi+2)) ...
          *(RIPP.FS*(1-RIPP.theta)*exp(-(v/RIPP.vs)^RIPP.ds)+RIPP.Fv*v)*u^(RIPP.chi+1) ;
end
FIwan = RIPP.F0+sign(v)*(RIPP.FS-sign(v)*RIPP.F0)/RIPP.FS*FIwan ;
RIPP.Flast = FIwan ;

force = Fpin + FIwan ;
