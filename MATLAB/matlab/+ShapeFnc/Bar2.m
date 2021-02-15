classdef Bar2 < ShapeFnc.Shape1D
    properties
        PltConn = [1 2]
    end
    
    methods
        %% All Interp matrices are nPoints x nNodes
        function obj=Bar2(conn,nodexyz)
            obj=obj@ShapeFnc.Shape1D(conn,nodexyz);
            obj.Dimension=1;
            obj.NumSides=1;
            obj.NumIntPts=2;
            obj.NumMassIntPts=obj.NumIntPts;
            obj.NumSurfIntPts=2;
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
        %%%
        function [conn,coord]=mkHigherOrder(obj)
            conn=obj.Connectivity;
            nelem=size(conn,1);
            coord=obj.Coordinates;
            n=size(obj.Coordinates,1);
            coord=[coord
                (coord(conn(:,1),:)+coord(conn(:,2),:))/2];
            conn=[conn ((n+1):(n+nelem))'];
        end
    end
    methods (Static)
        function N=Interpolate(pts)
            N=0.5*[(1-pts(:,1)) (1+pts(:,1))];
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
            dNx=0.5*[-ones(n,1) ones(n,1)];
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
            nodes=[1 2];
            dirs=2;
            node=nodes(sidenum);
            dir=dirs(sidenum)';
        end
        %%%%%%%%
        function [v,conn,elemtype]=GenSubCells(npts,type)
            %
            if nargin<2,
                type='s';
            end
            elemtype='bar';
            if strcmpi(type,'g'),
                if npts==1,
                        v=[-1;1];
                        conn=[1 2];
                elseif npts==2,
                    v=[-1;0;1];
                    conn=[1 2;2 3];
                else
                    [gpts]=ShapeFnc.Shape3D.GaussWeights(npts);
                    conn=[(1:npts)' (2:npts+1)'];
                    v=[-1;sum(gpts(conn(1:end-1,:)),2)/2;1];
                end
            else
               if npts==1,
                   v=[-1;1];
                   conn=[1 2];
               elseif npts==2,
                   v=[-1;0;1];
                   conn=[1 2;2 3];
               else
                   v=[-1;((1:(npts-1))*2/npts-1)';1];
                   conn=[(1:npts)' (2:npts+1)'];
               end
            end     
            %
        end
        %%
        function conn=ReverseElement(conn)
            conn=conn(:,[2 1]);
        end
        function conn=RandomizeConn(conn)
            m1=[1 2
                2 1];
            for i=1:size(conn,1),
                vi1=m1(randi(2,1),:);
                
                conn(i,:)=conn(i,vi1);
            end
        end
        %%%
    end
end
