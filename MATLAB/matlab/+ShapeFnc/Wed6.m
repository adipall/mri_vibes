classdef Wed6 < ShapeFnc.Shape3D
    properties
        PltConn = [1 2 3 1 4 5 6 4 6 3 2 5]
    end
    
    methods
        %% All Interp matrices are nPoints x nNodes
        function obj=Wed6(conn,nodexyz)
            obj=obj@ShapeFnc.Shape3D(conn,nodexyz);
            obj.NumSides=5;
            obj.NumIntPts=2;
            obj.NumMassIntPts=obj.NumIntPts;
            obj.NumSurfIntPts=2;
        end
        
        function [pts,weights,idnum]=getIntPts(obj,n)
            if nargin<2,
                n=obj.NumIntPts;
            end
            [pts,weights,idnum]=ShapeFnc.Shape3D.WedIntegration(n);
        end
        function [pts,weights]=getSurfIntPts(obj,side,n)
            if nargin==2,
                n=obj.NumSurfIntPts;
            end
            switch side
                case {1,2,3}
                    %% retangular surface
                    [pz,weights]=ShapeFnc.Shape3D.GaussWeights(n);
                    pt=(pz+1)/2; % scale between 0 and 1
                    wt=weights/2;
                    [pz,pt]=meshgrid(pz,pt);
                    [wts1,wts2]=meshgrid(weights,wt);
                    weights=wts1(:).*wts2(:);
                    switch side
                        case 1
                            pts=[pt(:) 1-pt(:) pz(:)];
                        case 2
                            pts=[zeros(numel(pt),1) pt(:) pz(:)];
                        case 3
                            pts=[1-pt(:) zeros(numel(pt),1) pz(:)];
                    end
                case {4,5}
                    %% triangular surface
                    [pt,weights]=ShapeFnc.Shape3D.TriIntegration(n);
                    switch side
                        case 4
                            pts=[pt -ones(size(pt,1),1)];
                        case 5
                            pts=[pt ones(size(pt,1),1)];
                    end
                otherwise
                    error('Unknown side in Wed6')
            end
        end
        function [pts,weights]=getGenIntPts(n,limits)
            if nargin==0,
                weights=ones(1,6)/6;
                a=0.5773502690;
                b=1/6;
                c=2/3;
                pts=[b b -a
                    c b -a
                    b c -a
                    b b a
                    c b a
                    b c a]';
            else
                if length(n)==1,
                    [pt,wt]=Shape3D.TriIntegration([n n]);
                    [pz,wz]=Shape3D.GaussWeights(n);
                elseif length(n)==3,
                    [pt,wt]=Shape3D.TriIntegration([n(1) n(2)]);
                    [pz,wz]=Shape3D.GaussWeights(n(3));
                else
                    error('Incorrect number of gauss points specified')
                end
            end
        end
        function [conn,coord]=mkHigherOrder(obj)
            conn=obj.Connectivity;
            nelem=size(conn,1);
            coord=obj.Coordinates;
            n=size(obj.Coordinates,1);
            coord=[coord
                (coord(conn(:,1),:)+coord(conn(:,2),:))/2
                (coord(conn(:,2),:)+coord(conn(:,3),:))/2
                (coord(conn(:,3),:)+coord(conn(:,1),:))/2
                (coord(conn(:,1),:)+coord(conn(:,4),:))/2
                (coord(conn(:,2),:)+coord(conn(:,5),:))/2
                (coord(conn(:,3),:)+coord(conn(:,6),:))/2
                (coord(conn(:,4),:)+coord(conn(:,5),:))/2
                (coord(conn(:,5),:)+coord(conn(:,6),:))/2
                (coord(conn(:,6),:)+coord(conn(:,4),:))/2];
            conn=[conn ((n+1):(n+nelem))' ((n+nelem+1):(n+2*nelem))' ...
                ((n+2*nelem+1):(n+3*nelem))' ((n+3*nelem+1):(n+4*nelem))' ...
                ((n+4*nelem+1):(n+5*nelem))' ((n+5*nelem+1):(n+6*nelem))' ...
                ((n+6*nelem+1):(n+7*nelem))' ((n+7*nelem+1):(n+8*nelem))' ...
                ((n+8*nelem+1):(n+9*nelem))'];
        end
    end
    methods (Static)
        function N=Interpolate(pts)
            % corner nodes
            a1=(1-pts(:,3));
            b1=(1+pts(:,3));
            t=1-pts(:,1)-pts(:,2);
            N=0.5*[a1.*pts(:,1) ...
                a1.*pts(:,2) ...
                a1.*t ...
                b1.*pts(:,1) ...
                b1.*pts(:,2) ...
                b1.*t];
        end
        function [pts]=LocalCentroid
            pts=[1/3 1/3 0];
        end
        function pts=LocalNode
          pts=[1  0   -1
              0   1   -1
              0   0   -1
              1   0    1
              0   1    1
              0   0    1];
        end
              
        function [xb,yb,zb]=LocalBounds
            xb=[0 1];
            yb=[0 1];
            zb=[-1 1];
        end
        %%%
        function dNx=lderivX(pts)
            a1=(1-pts(:,3));
            b1=(1+pts(:,3));
            z=zeros(size(pts,1),1);
            dNx=0.5*[a1 ...
                z ...
                -a1 ...
                b1 ...
                z ...
                -b1];
        end
        function dNy=lderivY(pts)
            a1=(1-pts(:,3));
            b1=(1+pts(:,3));
            z=zeros(size(pts,1),1);
            dNy=0.5*[z ...
                a1 ...
                -a1 ...
                z ...
                b1 ...
                -b1];
        end
        function dNz=lderivZ(pts)
            t=1-pts(:,1)-pts(:,2);
            dNz=0.5*[-pts(:,1) ...
                -pts(:,2) ...
                -t ...
                pts(:,1) ...
                pts(:,2) ...
                t];
        end
        %%%
        function [node,dir]=SideDef(sidenum)
            switch sidenum,
                case 1
                    node=[1,2,5,4];
                case 2
                    node=[2,3,6,5];
                case 3
                    node=[1,4,6,3];
                case 4
                    node=[1,3,2];
                case 5
                    node=[4,5,6]; % 3 node side
            end
            dirs=[4 1 2 -3 3];
            dir=dirs(sidenum)';
        end
        function [v,conn,elemtype]=GenSubCells(npts,type)
            if strcmpi(type,'g')
                if npts>5,
                    error('Wed6.GenSubCells not defined for gauss points greater than 5 order')
                end
            end
            [v,conn,elemtype]=ShapeFnc.Wed6.dbl_mesh(npts,type);
        end
        function [v,conn,elemtype]=dbl_mesh(npts,type)
            [vt,connt,elemtype]=ShapeFnc.Tri3.GenSubCells(npts,type);
            switch npts
                case 1,
                    npts=1;
                case 2
                    npts=2;
                case 3
                    npts=3;
                case 4
                    npts=3;
                case 5
                    npts=3;
                otherwise
                    error('number of npts not know')
            end
            [vz,connz]=ShapeFnc.Bar2.GenSubCells(npts,type);
            z=repmat(vz',size(vt,1),1);
            v=[repmat(vt,length(vz),1) z(:)];
            switch npts
                case 1
                    conn=[connt connt+size(vt,1)];
                case 2
                    conn=[connt connt+size(vt,1)
                        connt+size(vt,1) connt+2*size(vt,1)];
                case 3
                    conn=[connt connt+size(vt,1)
                        connt+size(vt,1) connt+2*size(vt,1)
                        connt+2*size(vt,1) connt+3*size(vt,1)];
                case 4
                    conn=[connt connt+size(vt,1)
                        connt+size(vt,1) connt+2*size(vt,1)
                        connt+2*size(vt,1) connt+3*size(vt,1)
                        connt+3*size(vt,1) connt+4*size(vt,1)];
                otherwise
                    error('Unknown number of points')
            end
        end
        function conn=ReverseElement(conn)
            conn=conn(:,[1 3 2 4 6 5]);
        end
        function conn=RandomizeConn(conn)
            m1=[1 2 3 4 5 6
                2 3 1 5 6 4
                3 1 2 6 4 5
                4 6 5 1 3 2
                6 5 4 3 2 1
                5 4 6 2 1 3];
            for i=1:size(conn,1),
                vi1=m1(randi(6,1),:);
                
                conn(i,:)=conn(i,vi1);
            end
        end
    end
end
