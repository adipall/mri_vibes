classdef Qua4 < ShapeFnc.Shape2D
    properties
        PltConn = [1 2 3 4 1]
    end
    
    methods
        %% All Interp matrices are nPoints x nNodes
        function obj=Qua4(conn,nodexyz)
            obj=obj@ShapeFnc.Shape2D(conn,nodexyz);
            obj.NumSides=4;
            obj.NumIntPts=2;
            obj.NumMassIntPts=obj.NumIntPts;
            obj.NumSurfIntPts=2;
        end
        function [pts,weights,gidx]=getIntPts(obj,n)
            if nargin==1,
                n=obj.NumSurfIntPts;
            end
            if length(n)==1,
                [pts,weights,gidx]=ShapeFnc.Shape3D.GaussWeights([n n]);
            elseif length(n)==2,
                [pts,weights,gidx]=ShapeFnc.Shape3D.GaussWeights(n);
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
        function [pts,weights]=getSurfIntPts(obj,side,n)
            [node,dir]=obj.SideDef(side);
            if nargin==2,
                n=obj.NumIntPts;
            end
            if length(n)==1,
                [pt,weights]=ShapeFnc.Shape3D.GaussWeights(n);
            else
                error('Incorrect number of gauss points specified')
            end
            switch abs(dir)
                case 1
                    pts=[sign(dir)*ones(size(pt,1),1) pt(:,1)];
                case 2
                    pts=[pt(:,1) sign(dir)*ones(size(pt,1),1)];
                otherwise
                    error('Unknown direction in Shape2D.')
            end
        end
        %%%
        function [PTSI]=GlobalIntercept(obj,PTS,vel,elemnum)
            if nargin==3,
                elemnum=1:size(obj.Connectivity,1);
            end
            mag=sum(vel.*vel,2);
            t=vel./mag(:,[1 1]);
            c=obj.Coordinates(:,1:2);
            conn=obj.Connectivity(elemnum,:);
            rc=[(c(conn(:,2),:)+c(conn(:,1),:)) ...
                (c(conn(:,3),:)+c(conn(:,2),:)) ...
                (c(conn(:,4),:)+c(conn(:,3),:)) ...
                (c(conn(:,1),:)+c(conn(:,4),:))]/2;
            n=[c(conn(:,2),:)-c(conn(:,1),:) ...
                c(conn(:,3),:)-c(conn(:,2),:) ...
                c(conn(:,4),:)-c(conn(:,3),:) ...
                c(conn(:,1),:)-c(conn(:,4),:)];
            nn=[sqrt(sum(n(:,[1 2]).*n(:,[1 2]),2)) ...
                sqrt(sum(n(:,[3 4]).*n(:,[3 4]),2)) ...
                sqrt(sum(n(:,[5 6]).*n(:,[5 6]),2)) ...
                sqrt(sum(n(:,[7 8]).*n(:,[7 8]),2))];
            n=[n(:,2)./nn(:,1) -n(:,1)./nn(:,1) ...
                n(:,4)./nn(:,2) -n(:,3)./nn(:,2) ...
                n(:,6)./nn(:,3) -n(:,5)./nn(:,3) ...
                n(:,8)./nn(:,4) -n(:,7)./nn(:,4)];

            r=rc-[PTS PTS PTS PTS];
            r=r.*n;
            num=[r(:,1)+r(:,2) r(:,3)+r(:,4) r(:,5)+r(:,6) r(:,7)+r(:,8)];
            t=n.*[t t t t];
            den=[t(:,1)+t(:,2) t(:,3)+t(:,4) t(:,5)+t(:,6) t(:,7)+t(:,8)];
            
            alpha=1./sort(den./num,2);
            switch 1
                case 1
                    %%for particle tracing
                    vi=4;
                case 2
                    %% for contact
                    vi=1;
            end
                    
            PTSI=PTS+[alpha(:,vi) alpha(:,vi)].*vel./mag(:,[1 1]);
           
        end
    end

    methods (Static)
        function N=Interpolate(pts)
            % corner nodes
            a1=(1-pts(:,1));a2=(1-pts(:,2));
            b1=(1+pts(:,1));b2=(1+pts(:,2));
            N=0.25*[a1 .* a2 ...
                b1 .* a2 ...
                b1 .* b2 ...
                a1 .* b2];
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
            a2=(1-pts(:,2));
            b2=(1+pts(:,2));
            dNx=0.25*[ -a2 ...
                a2 ...
                b2 ...
                -b2];
        end
        function dNy=lderivY(pts)
            a1=(1-pts(:,1));
            b1=(1+pts(:,1));
            dNy=0.25*[ -a1 ...
                -b1 ...
                b1 ...
                a1];
        end
        function dNz=lderivZ(pts)
            dNz=zeros(size(pts,1),4);
        end
        function [node,dir]=SideDef(sidenum)  
            nodes=[1,2
                2,3
                3,4
                4,1];
            dirs=[-2 1 2 -1];
            node=nodes(sidenum,:);
            dir=dirs(sidenum)';        
        end
        %%%%%%%
        function [v,conn,elemtype]=GenSubCells(npts,type)
            if nargin<2,
                type='s';
            end
            elemtype='quad';
            vx=ShapeFnc.Bar2.GenSubCells(npts,type);
            vy=vx;

            [vx,vy]=meshgrid(vx,vy);
            v=[vx(:) vy(:)];
            conn=ShapeFnc.Qua4.gen_conn2(length(vx),length(vy));
        end   
        function conn=gen_conn2(nx,ny)
            %% nx is the number of nodes in the x direction
            %% ny is the number of nodes in the y direction
            conn=zeros((nx-1)*(ny-1),4);
            for i=1:(nx-1),
                idx=((i-1)*(ny-1)+1):(i*(ny-1));
                startnode=((i-1)*ny+1:(i*ny-1))';
                conn(idx,:)=[startnode startnode+ny startnode+ny+1 startnode+1];
            end
        end
        %%%%%%%
        function conn=ReverseElement(conn)
            conn=conn(:,[1 4 3 2]);
        end
        function conn=RandomizeConn(conn)
            m1=[1 2 3 4
                2 3 4 1
                3 4 1 2
                4 1 2 3];
            for i=1:size(conn,1),
                vi1=m1(randi(4,1),:);
                
                conn(i,:)=conn(i,vi1);
            end
        end
    end
end
