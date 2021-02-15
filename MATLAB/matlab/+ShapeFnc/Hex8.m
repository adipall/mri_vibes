classdef Hex8 < ShapeFnc.Shape3D
    properties 
        PltConn = [1 2 6 5 1 4 3 7 8 4 8 5 6 7 3 2]
    end
    
    methods
        %% All Interp matrices are nPoints x nNodes
        function obj=Hex8(conn,nodexyz)
            obj=obj@ShapeFnc.Shape3D(conn,nodexyz);
            obj.NumSides=6;
            obj.NumIntPts=2;
            obj.NumMassIntPts=obj.NumIntPts;
            obj.NumSurfIntPts=2;
        end
        function [pts,weights,idnum]=getIntPts(obj,n)
            if nargin==1,
                n=obj.NumIntPts;
            end
            if length(n)==1,
                [pts,weights,idnum]=ShapeFnc.Shape3D.GaussWeights([n n n]);
            elseif length(n)==3,
                [pts,weights,idnum]=ShapeFnc.Shape3D.GaussWeights(n);
            else
                error('Incorrect number of gauss points specified')
            end
            gp_id=ShapeFnc.Shape3D.genGaussPointNames(pts);
        end
        function [pts,weights]=getGenIntPts(obj,n,limits)
            if nargin==1,
                n=obj.NumIntPts;
            end
            [pts,weights]=ShapeFnc.Shape3D.GeneralGaussWeights([n n n],limits);
        end
        function [pts,weights]=getSurfIntPts(obj,side,n)
            [node,dir]=obj.SideDef(side);
            if nargin==2,
                n=obj.NumIntPts;
            end
            if length(n)==1,
                [pt,weights]=ShapeFnc.Shape3D.GaussWeights([n n]);
            else
                error('Incorrect number of gauss points specified')
            end
            switch abs(dir)
                case 1
                    pts=[sign(dir)*ones(size(pt,1),1) pt];
                case 2
                    pts=[pt(:,1) sign(dir)*ones(size(pt,1),1) pt(:,2)];
                case 3
                    pts=[pt sign(dir)*ones(size(pt,1),1)];
                otherwise
                    error('Unknown direction in Shape3D.')
            end
        end
        %%%
        function [PTSI]=GlobalIntercept(obj,PTS,vel,elemnum)
            if nargin==3,
                elemnum=1:size(obj.Connectivity,1);
            end
            mag=sum(vel.*vel,2);
            t=vel./mag(:,[1 1]);
            conn=obj.Connectivity(elemnum,:);
            x=obj.Coordinates(:,1);
            y=obj.Coordinates(:,2);
            z=obj.Coordinates(:,3);
            
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
        function [conn,coord]=mkHigherOrder(obj)
            conn=obj.Connectivity;
            nelem=size(conn,1);
            coord=obj.Coordinates;
            n=size(obj.Coordinates,1);
            coord=[coord
                (coord(conn(:,1),:)+coord(conn(:,2),:))/2
                (coord(conn(:,2),:)+coord(conn(:,3),:))/2
                (coord(conn(:,3),:)+coord(conn(:,4),:))/2
                (coord(conn(:,4),:)+coord(conn(:,1),:))/2
                (coord(conn(:,1),:)+coord(conn(:,5),:))/2
                (coord(conn(:,2),:)+coord(conn(:,6),:))/2
                (coord(conn(:,3),:)+coord(conn(:,7),:))/2
                (coord(conn(:,4),:)+coord(conn(:,8),:))/2
                (coord(conn(:,5),:)+coord(conn(:,6),:))/2
                (coord(conn(:,6),:)+coord(conn(:,7),:))/2
                (coord(conn(:,7),:)+coord(conn(:,8),:))/2
                (coord(conn(:,8),:)+coord(conn(:,5),:))/2];
            conn=[conn ((n+1):(n+nelem))' ((n+nelem+1):(n+2*nelem))' ...
                ((n+2*nelem+1):(n+3*nelem))' ((n+3*nelem+1):(n+4*nelem))' ...
                ((n+4*nelem+1):(n+5*nelem))' ((n+5*nelem+1):(n+6*nelem))' ...
                ((n+6*nelem+1):(n+7*nelem))' ((n+7*nelem+1):(n+8*nelem))' ...
                ((n+8*nelem+1):(n+9*nelem))' ((n+9*nelem+1):(n+10*nelem))' ...
                ((n+10*nelem+1):(n+11*nelem))' ((n+11*nelem+1):(n+12*nelem))'];
        end
 
    end
    methods (Static)
        function N=Interpolate(pts)
            % corner nodes
            a1=(1-pts(:,1));a2=(1-pts(:,2));a3=(1-pts(:,3));
            b1=(1+pts(:,1));b2=(1+pts(:,2));b3=(1+pts(:,3));
            N=0.125*[a1 .* a2 .* a3 ...
                b1 .* a2 .* a3 ...
                b1 .* b2 .* a3 ...
                a1 .* b2 .* a3 ...
                a1 .* a2 .* b3 ...
                b1 .* a2 .* b3 ...
                b1 .* b2 .* b3 ...
                a1 .* b2 .* b3];
        end
        function [pts]=LocalCentroid
            pts=[0 0 0];
        end
        function pts=LocalNode
            pts=[-1 -1 -1
                1 -1 -1
                1  1 -1
                -1 1 -1
                -1 -1 1
                1 -1 1
                1 1 1
                -1 1 1];
        end
        function [xb,yb,zb]=LocalBounds
            xb=[-1 1];
            yb=[-1 1];
            zb=[-1 1];
        end
        
        function dNx=lderivX(pts)
            a2=(1-pts(:,2));a3=(1-pts(:,3));
            b2=(1+pts(:,2));b3=(1+pts(:,3));
            dNx=0.125*[-a2 .* a3, ...
                a2 .* a3, ...
                b2 .* a3, ...
                -b2 .* a3 ...
                -a2 .* b3, ...
                a2 .* b3, ...
                b2 .* b3, ...
                -b2 .* b3];
        end
        function dNy=lderivY(pts)
            a1=(1-pts(:,1));a3=(1-pts(:,3));
            b1=(1+pts(:,1));b3=(1+pts(:,3));
            dNy=0.125*[-a1 .* a3, ...
                -b1 .* a3, ...
                b1 .* a3, ...
                a1 .* a3 ...
                -a1 .* b3, ...
                -b1 .* b3, ...
                b1 .* b3, ...
                a1 .* b3];
        end
        function dNz=lderivZ(pts)
            a1=(1-pts(:,1));a2=(1-pts(:,2));
            b1=(1+pts(:,1));b2=(1+pts(:,2));
            dNz=0.125*[ -a1 .* a2, ...
                -b1 .* a2, ...
                -b1 .* b2, ...
                -a1 .* b2 ...
                a1 .* a2, ...
                b1 .* a2, ...
                b1 .* b2, ...
                a1 .* b2];
        end
        function [node,dir]=SideDef(sidenum)
            nodes=[1 2 6 5
                2 3 7 6
                3 4 8 7
                1 5 8 4
                1 4 3 2
                5 6 7 8];
            dirs=[-2 1 2 -1 -3 3];
            node=nodes(sidenum,:);
            dir=dirs(sidenum)';
        end
        %%%%%%%%
        function [v,conn,elemtype]=GenSubCells(npts,type)
            %
            if nargin<2,
                type='s';
            end
            elemtype='hex';
            %
            vx=ShapeFnc.Bar2.GenSubCells(npts,type);
            vy=vx;
            vz=vx;

            [vx,vy,vz]=meshgrid(vx,vy,vz);
            v=[vx(:),vy(:),vz(:)];
            conn=ShapeFnc.Hex8.gen_conn3(length(vx),length(vy),length(vz));
            
            
        end
        function conn=gen_conn3(nx,ny,nz)
            %% nx is the number of nodes in the x direction
            %% ny is the number of nodes in the y direction
            %% nz is the number of nodes in the z direction
            conn=zeros((nx-1)*(ny-1)*(nz-1),8);
            conn2=ShapeFnc.Qua4.gen_conn2(nx,ny);
            for i=1:(nz-1),
                idx=((i-1)*(nx-1)*(nz-1)+1):(i*(nx-1)*(ny-1));
                conn(idx,:)=[(i-1)*nx*ny+conn2 conn2+i*nx*ny];
            end
        end
        %%%%%%%%
        function conn=ReverseElement(conn)
            conn=conn(:,[1 4 3 2 5 8 7 6]);
        end
        function conn=RandomizeConn(conn)
            m1=[1 2 3 4
                2 3 4 1
                3 4 1 2
                4 1 2 3];
            m2=[1 2];
            for i=1:size(conn,1),
                vi1=m1(randi(4,1),:);
                vi2=m2(randi(2,1));
                vi3=randi(3,1);
                
                switch vi3
                    case 1
                        vi=1:8;
                    case 2
                        vi=[3 2 6 7 4 1 5 8];
                    case 3
                        vi=[1 5 6 2 4 8 7 3];
                end
                conn(i,:)=conn(i,vi);
                if vi2==2,
                    vi1=vi1(4:-1:1);
                    vi=[vi1+4 vi1];
                else
                    vi=[vi1 vi1+4];
                end
                conn(i,:)=conn(i,vi);
            end
        end
    end
end
