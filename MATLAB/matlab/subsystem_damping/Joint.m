classdef Joint
    % This is an abstract class that is never called
    %
    properties % these values must be updated by updateForce and updateStiffness
        Kt  % current value of tangent stiffness
        K0  % no load stiffness
        Force % current value of the force in the iwan element
    end
    
    methods
        function obj=Joint(params)
            if nargin>0,
                obj.Kt=params;
            end
        end
    end
    methods (Abstract)
        [obj,force]=updateForce(obj,x,v)
        [obj,stiff]=updateStiffness(obj,x,v)
        [stiff]=calcStiffness(obj,x,v)
        [obj]=updateStates(obj)
    end
end
