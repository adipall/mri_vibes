classdef GlobalVars < handle
    properties
        Name
        Data
    end
    methods
        function obj=GlobalVars(name)
            if nargin>0,
                for i=1:length(name),
                    obj(i)=FEMesh.GlobalVars;
                end
                for i=1:length(name),
                    obj(i).Name=name{i};
                end
            end
        end
        function flag=isGlobalVar(obj,gname)
            flag=0;
            for j=1:size(obj,2),
                if strcmpi(obj(1,j).Name,gname)
                    flag=j;
                    break;
                end
            end
        end
    end
end
