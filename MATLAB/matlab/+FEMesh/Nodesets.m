classdef Nodesets
    properties
        ID
        Nodes
        DistFactors
        Status  = 1
        Name
    end
    methods
        function obj=Nodesets(id)
            if nargin>0,
                for i=1:length(id),
                    obj(i)=FEMesh.Nodesets;
                end
                for i=1:length(id),
                    obj(i).ID=id(i);
                end
            end
        end
        function [idx]=id2idx(obj,id)
            ids=[obj.ID];
            idx=find(ids==id);
        end
        function [ids]=getNodesetIDs(obj)
            ids=[obj.ID];
        end
    end
end