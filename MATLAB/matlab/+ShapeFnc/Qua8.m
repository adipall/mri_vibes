classdef Qua8 < ShapeFnc.Shape2D
    properties
        PltConn = [1 5 2 6 3 7 4 8 1]
    end
    
    methods
        %% All Interp matrices are nPoints x nNodes
        function obj=Qua8(conn,nodexyz)
            obj=obj@ShapeFnc.Shape2D(conn,nodexyz);
            obj.NumSides=4;
            obj.NumIntPts=3;
            obj.NumMassIntPts=obj.NumIntPts;
            obj.NumSurfIntPts=3;
        end
        function [pts,weights]=getIntPts(obj,n)
            if nargin==1,
                n=obj.NumIntPts;
            end
            if length(n)==1,
                [pts,weights]=ShapeFnc.Shape3D.GaussWeights([n n]);
            elseif length(n)==2,
                [pts,weights]=ShapeFnc.Shape3D.GaussWeights(n);
            else
                error('Incorrect number of gauss points specified')
            end
        end
        function [pts,weights]=getGenIntPts(obj,n,limits)
            if nargin==1,
                n=obj.NumIntPts;
            end
            [pts,weights]=ShapeFnc.Shape3D.GeneralGaussWeights([n n],limits);
        end
        function [pts,weights]=getSurfIntPts(obj,n)
            [node,dir]=obj.SideDef(side);
            if nargin==1,
                n=obj.NumIntPts;
            end
            if length(n)==1,
                [pt,weights]=ShapeFnc.Shape3D.GaussWeights(n);
            else
                error('Incorrect number of gauss points specified')
            end
            switch abs(dir)
                case 1
                    pts=[sign(dir)*ones(size(pt,1),1) pt];
                case 2
                    pts=[pt sign(dir)*ones(size(pt,1),1)];
                otherwise
                    error('Unknown direction in Shape2D.')
            end
        end
    end
    methods(Static)
        function N=Interpolate(pts)
            % corner nodes
            N=zeros(size(pts,1),8);
            a1=(1-pts(:,1));a2=(1-pts(:,2));
            b1=(1+pts(:,1));b2=(1+pts(:,2));
            N(:,1:4)=0.25*[a1 .* a2 .* (-pts(:,1) - pts(:,2) - 1) ...
                b1 .* a2 .* (pts(:,1) - pts(:,2) - 1) ...
                b1 .* b2 .* (pts(:,1) + pts(:,2) - 1) ...
                a1 .* b2 .* (-pts(:,1) + pts(:,2) - 1)];
            %
            d=(1-pts(:,1).*pts(:,1));
            N(:,5)=0.5*d .*a2;
            N(:,7)=0.5*d.*b2;
            %
            d=(1-pts(:,2).*pts(:,2));
            N(:,6)=0.5*b1.*d;
            N(:,8)=0.5*a1.*d;
        end
        function [pts]=LocalCentroid
            pts=[0 0];
        end
        function [xb,yb,zb]=LocalBounds
            xb=[-1 1];
            yb=[-1 1];
            zb=[0 0];
        end
        
        function dNx=lderivX(pts)
            dNx=zeros(size(pts,1),8);
            a1=(1-pts(:,1));a2=(1-pts(:,2));
            b1=(1+pts(:,1));b2=(1+pts(:,2));
            dNx(:,1:4)=0.25*[ -a2 .* (-pts(:,1) - pts(:,2) - 1) - a1 .* a2...
                a2 .* (pts(:,1) - pts(:,2) - 1) + b1 .* a2...
                b2 .* (pts(:,1) + pts(:,2) - 1) + b1 .* b2...
                -b2.* (-pts(:,1) + pts(:,2) - 1) - a1 .* b2 ];
            d=-2*pts(:,1);
            dNx(:,5)=0.5*d .*a2;
            dNx(:,7)=0.5*d.*b2;
            %
            d=(1-pts(:,2).*pts(:,2));
            dNx(:,6)=0.5*d;
            dNx(:,8)=-0.5*d;
        end
        function dNy=lderivY(pts)
            dNy=zeros(size(pts,1),8);
            a1=(1-pts(:,1));a2=(1-pts(:,2));
            b1=(1+pts(:,1));b2=(1+pts(:,2));
            dNy(:,1:4)=0.25*[ -a1 .* (-pts(:,1) - pts(:,2) - 1) - a1 .* a2...
                -b1 .* (pts(:,1) - pts(:,2) - 1) - b1 .* a2...
                b1 .* (pts(:,1) + pts(:,2) - 1) + b1 .* b2...
                a1.* (-pts(:,1) + pts(:,2) - 1)+ a1 .* b2];
            %
            d=(1-pts(:,1).*pts(:,1));
            dNy(:,5)=-0.5*d;
            dNy(:,7)=0.5*d;
            %
            d=-2*pts(:,2);
            dNy(:,6)=0.5*b1.*d;
            dNy(:,8)=0.5*a1.*d;
        end
        function dNz=lderivZ(pts)
            dNz=zeros(size(pts,1),8);
        end
        function [node,dir]=SideDef(sidenum)
            nodes=[1,2,5
                2,3,6
                3,4,7
                4,1,8];
            dirs=[-2 1 2 -1];
            node=nodes(sidenum,:);
            dir=dirs(sidenum)';
        end
        function [v,conn,elemtype]=GenSubCells(npts,type)
            if nargin<2,
                type='s';
            end
            tol=1e-14;
            [v,conn,elemtype]=ShapeFnc.Qua4.GenSubCells(npts,type);
            nv=size(v,1);
            v=[v
                (v(conn(:,2),:)+v(conn(:,1),:))/2
                (v(conn(:,3),:)+v(conn(:,2),:))/2
                (v(conn(:,4),:)+v(conn(:,3),:))/2
                (v(conn(:,4),:)+v(conn(:,1),:))/2];
            nc=size(conn,1);
            conn=[conn ...
                (nv+(1:nc))'...
                (nv+((nc+1):(2*nc)))' ...
                (nv+((2*nc+1):(3*nc)))' ...
                (nv+((3*nc+1):(4*nc)))'];
            size(v)
            [v,I,J]=unique(tol*round(v/tol),'rows');
            
            conn=J(conn);
            if npts==1,
                conn=conn';
            end
        end
        function conn=ReverseElement(conn)
            conn=conn(:,[1 4 3 2 8 7 6 5]);
        end
        function conn=RandomizeConn(conn)
            m1=[1 2 3 4 5 6 7 8
                2 3 4 1 6 7 8 5
                3 4 1 2 7 8 5 6
                4 1 2 3 8 5 6 7];
            for i=1:size(conn,1),
                vi1=m1(randi(4,1),:);
                
                conn(i,:)=conn(i,vi1);
            end
        end
    end
end
