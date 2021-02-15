classdef Nodes < handle
    properties
        Names
        Coordinates
        NodeNumMap
    end
    methods
        function obj=Nodes(coord,nmap,names)
            if nargin>0,
                obj.Coordinates=coord;
                obj.NodeNumMap=nmap(:);
                if nargin==3,
                    obj.Names=names;
                else
                    switch size(coord,2),
                        case 1
                            obj.Names={'x_coord'};
                        case 2
                            obj.Names={'x_coord','y_coord'};
                        case 3
                            obj.Names={'x_coord','y_coord','z_coord'};
                        otherwise
                            error('Unknown Dimension')
                    end
                end
            end
        end
        function display1(obj)
            
            [n,m]=size(obj.Coordinates);
            
            disp(sprintf(' NumNodes    NumDim'))
            disp(sprintf('%d          %d',n,m))
        end
    end
end
