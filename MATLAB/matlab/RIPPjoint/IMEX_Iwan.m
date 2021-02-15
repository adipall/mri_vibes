function [force, params] = IMEX_Iwan(u, params)

[force, ~, ~, joint_states_temp] ...
      = iwan_function_preferred(u,0,params.y0,params.joint);

% Update the joint_states vector
params.y0 = joint_states_temp ;
