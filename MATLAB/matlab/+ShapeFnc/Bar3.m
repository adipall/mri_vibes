classdef Bar3 < ShapeFnc.Shape3D
    properties
        PltConn = [1 3 2]
    end
    
    methods
        %% All Interp matrices are nPoints x nNodes
        function obj=Bar3(conn,nodexyz)
            obj=obj@ShapeFnc.Shape3D(conn,nodexyz);
            obj.Dimension=1;
            obj.NumSides=1;
            obj.NumIntPts=3;
            obj.NumMassIntPts=obj.NumIntPts;
            obj.NumSurfIntPts=3;
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
    end
    methods (Static)
        function N=Interpolate(pts)
            N=[0.5*(1-pts(:,1)).*pts(:,1) 0.5*(1+pts(:,1)).*pts(:,1) (1-pts(:,1).*pts(:,1))];
        end
        function [pts]=LocalCentroid
            pts=0;
        end
        function [xb,yb,zb]=LocalBounds
            xb=[-1 1];
            yb=[0 0];
            zb=[0 0];
        end
        
        function dNx=lderivX(pts)
            n=size(pts,1);
            dNx=[-0.5*ones(n,1) 0.5*ones(n,1) -2*ones(n,1)];
        end
        function dNy=lderivY(pts)
            dNy=zeros(size(pts,1),2);
        end
        function dNz=lderivZ(pts)
            dNz=zeros(size(pts,1),2);
        end
        function [node,dir]=SideDef(sidenum)
            if any(sidenum>1),
                error('Only one side defined for a bar2 element')
            end
            nodes=[1 2 3];
            dirs=2;
            node=nodes(sidenum);
            dir=dirs(sidenum)';
        end
        function [v,conn,elemtype]=GenSubCells(npts,type)
            if nargin<2,
                type='s';
            end
            tol=1e-14;
            [v,conn,elemtype]=ShapeFnc.Bar2.GenSubCells(npts,type);
            nv=size(v,1);
            v=[v
                (v(conn(:,2),:)+v(conn(:,1),:))/2];
            nc=size(conn,1);
            conn=[conn ...
                (nv+(1:nc))'];
            [v,I,J]=unique(tol*round(v/tol),'rows');
            conn=J(conn);
            if npts==1,
                conn=conn';
            end
        end
        %%
        function conn=ReverseElement(conn)
            conn=conn(:,[2 1 3]);
        end
        function conn=RandomizeConn(conn)
            m1=[1 2 3
                2 1 3];
            for i=1:size(conn,1),
                vi1=m1(randi(2,1),:);
                
                conn(i,:)=conn(i,vi1);
            end
        end
    end
end
