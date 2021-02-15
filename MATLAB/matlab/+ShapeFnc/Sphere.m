classdef Sphere < ShapeFnc.Shape1D
    properties
        PltConn = [1]
    end
    
    methods
        %% All Interp matrices are nPoints x nNodes
        function obj=Sphere(conn,nodexyz)
            obj=obj@ShapeFnc.Shape1D(conn,nodexyz);
            obj.Dimension=1;
            obj.NumSides=1;
            obj.NumIntPts=1;
            obj.NumMassIntPts=obj.NumIntPts;
            obj.NumSurfIntPts=1;
        end
        function [pts,weights,idnum]=getIntPts(obj,n)
            if nargin==1,
                n=obj.NumIntPts;
            end
            if length(n)==1,
                [pts,weights,idnum]=ShapeFnc.Shape3D.GaussWeights(n);
            else
                error('Incorrect number of gauss points specified')
            end
        end
        function [pts,weights]=getGenIntPts(obj,n,limits)
            if nargin==1,
                n=obj.NumIntPts;
            end
            [pts,weights]=ShapeFnc.Shape3D.GeneralGaussWeights(n,limits);
        end
        function [pts,weights]=getSurfIntPts(obj,n)
            if nargin==1,
                n=obj.NumSurfIntPts;
            end
            if length(n)==1,
                [pts,weights]=ShapeFnc.Shape3D.GaussWeights(n);
            else
                error('Incorrect number of gauss points specified')
            end
        end
        function [vol]=Volume(obj,elemnum)
            if nargin==1,
                elemnum=1:size(obj.Connectivity,1);
            end
            nelem=length(elemnum);
            vol=ones(nelem,1);
        end
    end
    methods (Static)
        function N=Interpolate(pts)
            N=ones(size(pts,1),1);
        end
        function [pts]=LocalCentroid
            pts=[0];
        end
        function [xb,yb,zb]=LocalBounds
            xb=[-1 1];
            yb=[0 0];
            zb=[0 0];
        end
        
        function dNx=lderivX(pts)
            n=size(pts,1);
            dNx=zeros(n,1);
        end
        function dNy=lderivY(pts)
            dNy=zeros(size(pts,1),2);
        end
        function dNz=lderivZ(pts)
            dNz=zeros(size(pts,1),2);
        end
        function [node,dir]=SideDef(sidenum)
            if any(sidenum>0),
                error('No sides are defined for a sphere element')
            end
            
        end
        %%%%%%%%
        function [v,conn,elemtype]=GenSubCells(npts,type)
            %
            error('Sphere elements cannot be subdivided')
            %
        end
        %%
       function conn=ReverseElement(conn)
            conn=conn;
        end
        function conn=RandomizeConn(conn)
            conn=conn;
        end
    end
end
