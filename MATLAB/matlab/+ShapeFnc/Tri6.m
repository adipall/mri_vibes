classdef Tri6 < ShapeFnc.Shape2D
    properties
        PltConn = [1 4 2 5 3 6 1]
    end
    
    methods
        %% All Interp matrices are nPoints x nNodes
        function obj=Tri6(conn,nodexyz)
            obj=obj@ShapeFnc.Shape2D(conn,nodexyz);
            obj.NumSides=3;
            obj.NumIntPts=3;
            obj.NumMassIntPts=obj.NumIntPts;
            obj.NumSurfIntPts=3;
        end
        function [pts,weights]=getIntPts(obj,n)
            if nargin==1,
                n=obj.NumIntPts;
            end
            if length(n)==1,
                [pts,weights]=ShapeFnc.Shape3D.TriIntegration(n);
            else
                error('Incorrect number of gauss points specified')
            end
        end
        function [pts,weights]=getGenIntPts(obj,n,limits)
            if nargin==1,
                n=3;
            end
            if length(n)==1,
                [pts,weights]=ShapeFnc.Shape3D.TriIntegration(n);
            else
                error('Incorrect number of gauss points specified')
            end
        end
        function [pts,weights]=getSurfIntPts(obj,side,n)
            if nargin==2,
                n=obj.NumSurfIntPts;
            end
            if length(n)==1,
                [pt,weights]=ShapeFnc.Shape3D.GaussWeights(n);
            else
                error('Incorrect number of gauss points specified')
            end
            pt=(pt+1)/2;  % scale between 0 and 1
            weights=weights/2;
            switch side
                case 1
                    pts=[pt 1-pt];
                case 2
                    pts=[zeros(size(pt,1),1) pt ];
                case 3
                    pts=[1-pt zeros(size(pt,1),1)];
                otherwise
                    error('Unknown side in Tri3')
            end
        end
    end
    methods (Static)
        function N=Interpolate(pts)
            % use area coordinates
            % corner nodes
            l3=1-pts(:,1)-pts(:,2);
            N=[(2*pts(:,1)-1).*pts(:,1) ...
                (2*pts(:,2)-1).*pts(:,2) ...
                (2*l3-1).*l3 ...
                4*pts(:,1).*pts(:,2) ...
                4*pts(:,2).*l3 ...
                4*pts(:,1).*l3];
        end
        function [pts]=LocalCentroid
            pts=[1/3 1/3];
        end
        function [xb,yb,zb]=LocalBounds
            xb=[0 1];
            yb=[0 1];
            zb=[0 1];
        end
        %%%
        function dNx=lderivX(pts)
            l3=1-pts(:,1)-pts(:,2);
            dNx=[4*pts(:,1)-1 ...
                zeros(size(pts,1),1) ...
                1-4*l3 ...
                4*pts(:,2) ...
                -4*pts(:,2) ...
                4*(l3-pts(:,1))];
        end
        function dNy=lderivY(pts)
            l3=1-pts(:,1)-pts(:,2);
            dNy=[zeros(size(pts,1),1) ...
                4*pts(:,2)-1 ...
                1-4*l3 ...
                4*pts(:,1) ...
                4*(l3-pts(:,2)) ...
                -4*pts(:,1)];
        end
        function dNz=lderivZ(pts)
            dNz=zeros(size(pts,1),6);
        end
        function [node,dir]=SideDef(sidenum)
            nodes=[1,2,4
                2,3,5
                3,1,6];
            dirs=[3 1 2];
            node=nodes(sidenum,:);
            dir=dirs(sidenum)';
        end
        function [v,conn]=GenSubCells(npts,type)
            if nargin<2,
                type='s';
            end
            tol=1e-14;
            [v,conn]=ShapeFnc.Tri3.GenSubCells(npts,type);
            nv=size(v,1);
            v=[v
                (v(conn(:,2),:)+v(conn(:,1),:))/2
                (v(conn(:,3),:)+v(conn(:,2),:))/2
                (v(conn(:,3),:)+v(conn(:,1),:))/2];
            nc=size(conn,1);
            conn=[conn ...
                (nv+(1:nc))'...
                (nv+((nc+1):(2*nc)))' ...
                (nv+((2*nc+1):(3*nc)))'];
            [v,I,J]=unique(tol*round(v/tol),'rows');
            conn=J(conn);
            if npts==1,
                conn=conn';
            end
        end
        function conn=ReverseElement(conn)
            conn=conn(:,[1 3 2 6 5 4]);
        end
        function conn=RandomizeConn(conn)
            m1=[1 2 3 4 5 6
                2 3 1 5 6 4
                3 1 2 6 4 5];
            for i=1:size(conn,1),
                vi1=m1(randi(3,1),:);
                
                conn(i,:)=conn(i,vi1);
            end
        end
    end
end
