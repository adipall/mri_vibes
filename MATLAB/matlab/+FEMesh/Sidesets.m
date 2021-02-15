classdef Sidesets
    properties
        ID
        Name
        Elements
        Sides
        DistFactors
        Status  = 1
        Nodes
        NumNodes
    end
   
    methods
        function obj=Sidesets(id)
            if nargin>0,
                for i=1:length(id)
                    obj(i)=FEMesh.Sidesets;
                end
                for i=1:length(id)
                    obj(i).ID=id(i);
                end
            end
        end
        function [idx]=id2idx(obj,id)
            ids=[obj.ID];
            idx=find(ids==id);
        end
        function [ids]=getSidesetIDs(obj)
            ids=[obj.ID];
        end
    end
 
end
