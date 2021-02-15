classdef LinDamp < Joint
    
    properties (Access=private)
        Zeta = 0  % % of critical damping
    end
    
    methods
        function obj=LinDamp(params)
            if nargin>0,
                % for modal damping in the Zeta=2*zeta*wn
                obj.Zeta=params(1);
                obj.Kt=params(2);
                obj=obj.updateStiffness(0,0);
                obj.K0=obj.Kt;
            end
        end
        function [obj,force]=updateForce(obj,x,v)
           obj.Force=obj.Zeta*v;
           force=obj.Force;
        end
        
        function [obj,kt]=updateStiffness(obj,x,v)
            % stiffness is independent of differential stiffness (x) and
            % Differential velocity (v)
            kt=obj.Kt;
        end
        function kt=calcStiffness(obj,x,v)
            kt=obj.Kt;
        end
        function obj=updateStates(obj)
            %% nothing to do here
        end
    end
end
