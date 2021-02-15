classdef NodalVars < handle
    properties
        Name
        Data
    end
    methods
        function obj=NodalVars(name)
            if nargin>0,
                for i=1:length(name),
                    obj(i)=FEMesh.NodalVars;
                end
                for i=1:length(name),
                    obj(i).Name=name{i};
                end
            end
        end
        function nvar=getnvar(obj,nname)
            for i=1:length(obj),
                if strcmp(nname,obj(i).Name)
                    nvar=obj(i).Data;
                    return
                end
            end
            warning('junk','Nodal Variable %s not found',nname)
            nvar=[];
        end
        function flag=isNodalVar(obj,nname)
            flag=0;
            for j=1:length(obj),
                if strcmpi(obj(j).Name,nname)
                    flag=j;
                    break;
                end
            end
        end
    end
end
