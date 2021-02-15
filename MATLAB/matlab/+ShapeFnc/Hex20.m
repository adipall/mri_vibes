classdef Hex20 < ShapeFnc.Shape3D
    properties 
        PltConn = [1 9 2 14 6 17 5 13 1 12 4 11 3 15 7 19 8 16 4 16 8 20 5 17 6 18 7 15 3 10 2]
    end

    methods
        %% All Interp matrices are nPoints x nNodes
        function obj=Hex20(conn,nodexyz)
            obj=obj@ShapeFnc.Shape3D(conn,nodexyz);
            obj.NumSides=6;
            obj.NumIntPts=3;
            obj.NumMassIntPts=obj.NumIntPts;
            obj.NumSurfIntPts=3;
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
        end
        function [pts,weights]=getGenIntPts(obj,n,limits)
            if nargin==1,
                n=obj.NumIntPts;
            end
            [pts,weights]=ShapeFnc.Shape3D.GeneralGaussWeights([n n n],limits);
        end
        function [pts,weights]=getSurfIntPts(obj,side,n)
            [node,dir]=obj.SideDef(side);
            if nargin==1,
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
    end
    methods (Static)
        function N=Interpolate(pts)
            N=zeros(size(pts,1),20);
            % corner nodes
            a1=(1-pts(:,1));a2=(1-pts(:,2));a3=(1-pts(:,3));
            b1=(1+pts(:,1));b2=(1+pts(:,2));b3=(1+pts(:,3));
            N(:,1:8)=0.125*[a1 .* a2 .* a3.*(-pts(:,1)-pts(:,2)-pts(:,3)-2) ...
                b1 .* a2 .* a3.*(pts(:,1)-pts(:,2)-pts(:,3)-2) ...
                b1 .* b2 .* a3.*(pts(:,1)+pts(:,2)-pts(:,3)-2) ...
                a1 .* b2 .* a3.*(-pts(:,1)+pts(:,2)-pts(:,3)-2) ...
                a1 .* a2 .* b3.*(-pts(:,1)-pts(:,2)+pts(:,3)-2) ...
                b1 .* a2 .* b3.*(pts(:,1)-pts(:,2)+pts(:,3)-2) ...
                b1 .* b2 .* b3.*(pts(:,1)+pts(:,2)+pts(:,3)-2) ...
                a1 .* b2 .* b3.*(-pts(:,1)+pts(:,2)+pts(:,3)-2)];
            d=(1-pts(:,1).*pts(:,1));
            N(:,9)=0.25*d.*a2.*a3;
            N(:,11)=0.25*d.*b2.*a3;
            N(:,17)=0.25*d.*a2.*b3;
            N(:,19)=0.25*d.*b2.*b3;
            %N(:,vi)=0.25*(1-repmat(x(:).*x(:),1,4)).*(1+b(:,vi)).*(1+c(:,vi));
            d=(1-pts(:,2).*pts(:,2));
            N(:,10)=0.25*b1.*d.*a3;
            N(:,12)=0.25*a1.*d.*a3;
            N(:,18)=0.25*b1.*d.*b3;
            N(:,20)=0.25*a1.*d.*b3;
            %N(:,vi)=0.25*(1+a(:,vi)).*(1-repmat(y(:).*y(:),1,4)).*(1+c(:,vi));
            d=(1-pts(:,3).*pts(:,3));
            N(:,13)=0.25*a1.*a2.*d;
            N(:,14)=0.25*b1.*a2.*d;
            N(:,15)=0.25*b1.*b2.*d;
            N(:,16)=0.25*a1.*b2.*d;
            %N(:,vi)=0.25*(1+a(:,vi)).*(1+b(:,vi)).*(1-repmat(z(:).*z(:),1,4));
            %
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
                -1 1 1
                0 -1 -1
                1 0 -1
                0 1 -1
                -1 0 -1
                -1 -1 0
                1 -1 0
                1 1 0
                -1 1 0
                0 -1 1
                1 0 1
                0 1 1
                -1 0 1];
        end
        function [xb,yb,zb]=LocalBounds
            xb=[-1 1];
            yb=[-1 1];
            zb=[-1 1];
        end

        function dNx=lderivX(pts)
            dNx=zeros(size(pts,1),20);
            a1=(1-pts(:,1));a2=(1-pts(:,2));a3=(1-pts(:,3));
            b1=(1+pts(:,1));b2=(1+pts(:,2));b3=(1+pts(:,3));
            dNx(:,1:8)=0.125*[-a2 .* a3 .*(-pts(:,1)-pts(:,2)-pts(:,3)-2) - a1 .* a2 .* a3 ...
                a2 .* a3 .*(pts(:,1)-pts(:,2)-pts(:,3)-2) + b1 .* a2 .* a3 ...
                b2 .* a3 .*(pts(:,1)+pts(:,2)-pts(:,3)-2) + b1 .* b2 .* a3 ...
                -b2 .* a3 .*(-pts(:,1)+pts(:,2)-pts(:,3)-2) - a1 .* b2 .* a3 ...
                -a2 .* b3 .*(-pts(:,1)-pts(:,2)+pts(:,3)-2) - a1 .* a2 .* b3 ...
                a2 .* b3 .*(pts(:,1)-pts(:,2)+pts(:,3)-2) + b1 .* a2 .* b3 ...
                b2 .* b3 .*(pts(:,1)+pts(:,2)+pts(:,3)-2) + b1 .* b2 .* b3 ...
                -b2 .* b3.*(-pts(:,1)+pts(:,2)+pts(:,3)-2) - a1 .* b2 .* b3];
            % Midside nodes
            d=-2*pts(:,1);
            dNx(:,9)=0.25*d.*a2.*a3;
            dNx(:,11)=0.25*d.*b2.*a3;
            dNx(:,17)=0.25*d.*a2.*b3;
            dNx(:,19)=0.25*d.*b2.*b3;
            %
            d=(1-pts(:,2).*pts(:,2));
            dNx(:,10)=0.25*d.*a3;
            dNx(:,12)=-0.25*d.*a3;
            dNx(:,18)=0.25*d.*b3;
            dNx(:,20)=-0.25*d.*b3;
            %
            d=(1-pts(:,3).*pts(:,3));
            dNx(:,13)=-0.25*a2.*d;
            dNx(:,14)=0.25*a2.*d;
            dNx(:,15)=0.25*b2.*d;
            dNx(:,16)=-0.25*b2.*d;
        end
        function dNy=lderivY(pts)
            dNy=zeros(size(pts,1),20);
            a1=(1-pts(:,1));a2=(1-pts(:,2));a3=(1-pts(:,3));
            b1=(1+pts(:,1));b2=(1+pts(:,2));b3=(1+pts(:,3));
            dNy(:,1:8)=0.125*[-a1 .* a3 .*(-pts(:,1)-pts(:,2)-pts(:,3)-2) - a1 .* a2 .* a3...
                -b1 .* a3.*(pts(:,1)-pts(:,2)-pts(:,3)-2) - b1 .* a2 .* a3...
                b1 .* a3 .*(pts(:,1)+pts(:,2)-pts(:,3)-2)+ b1 .* b2 .* a3...
                a1 .* a3 .*(-pts(:,1)+pts(:,2)-pts(:,3)-2) + a1 .* b2 .* a3...
                -a1 .* b3 .*(-pts(:,1)-pts(:,2)+pts(:,3)-2)- a1 .* a2 .* b3 ...
                -b1 .* b3 .*(pts(:,1)-pts(:,2)+pts(:,3)-2) - b1 .* a2 .* b3...
                b1 .* b3 .*(pts(:,1)+pts(:,2)+pts(:,3)-2) + b1 .* b2 .* b3...
                a1 .* b3.*(-pts(:,1)+pts(:,2)+pts(:,3)-2)+ a1 .* b2 .* b3];
            %
            d=(1-pts(:,1).*pts(:,1));
            dNy(:,9)=-0.25*d.*a3;
            dNy(:,11)=0.25*d.*a3;
            dNy(:,17)=-0.25*d.*b3;
            dNy(:,19)=0.25*d.*b3;
            %
            d=-2*pts(:,2);
            dNy(:,10)=0.25*b1.*d.*a3;
            dNy(:,12)=0.25*a1.*d.*a3;
            dNy(:,18)=0.25*b1.*d.*b3;
            dNy(:,20)=0.25*a1.*d.*b3;
            %
            d=(1-pts(:,3).*pts(:,3));
            dNy(:,13)=-0.25*a1.*d;
            dNy(:,14)=-0.25*b1.*d;
            dNy(:,15)=0.25*b1.*d;
            dNy(:,16)=0.25*a1.*d;
        end
        function dNz=lderivZ(pts)
            dNz=zeros(size(pts,1),20);
            a1=(1-pts(:,1));a2=(1-pts(:,2));a3=(1-pts(:,3));
            b1=(1+pts(:,1));b2=(1+pts(:,2));b3=(1+pts(:,3));
            dNz(:,1:8)=0.125*[ -a1 .* a2 .*(-pts(:,1)-pts(:,2)-pts(:,3)-2)- a1 .* a2 .* a3...
                -b1 .* a2 .*(pts(:,1)-pts(:,2)-pts(:,3)-2)-b1 .* a2 .* a3...
                -b1 .* b2 .*(pts(:,1)+pts(:,2)-pts(:,3)-2)- b1 .* b2 .* a3...
                -a1 .* b2 .*(-pts(:,1)+pts(:,2)-pts(:,3)-2)- a1 .* b2 .* a3...
                a1 .* a2  .*(-pts(:,1)-pts(:,2)+pts(:,3)-2)+ a1 .* a2 .* b3...
                b1 .* a2 .*(pts(:,1)-pts(:,2)+pts(:,3)-2) + b1 .* a2 .* b3...
                b1 .* b2 .*(pts(:,1)+pts(:,2)+pts(:,3)-2) + b1 .* b2 .* b3...
                a1 .* b2.*(-pts(:,1)+pts(:,2)+pts(:,3)-2)+ a1 .* b2 .* b3];
            %
            d=(1-pts(:,1).*pts(:,1));
            dNz(:,9)=-0.25*d.*a2;
            dNz(:,11)=-0.25*d.*b2;
            dNz(:,17)=0.25*d.*a2;
            dNz(:,19)=0.25*d.*b2;
            %
            d=(1-pts(:,2).*pts(:,2));
            dNz(:,10)=-0.25*b1.*d;
            dNz(:,12)=-0.25*a1.*d;
            dNz(:,18)=0.25*b1.*d;
            dNz(:,20)=0.25*a1.*d;
            %
            d=-2*pts(:,3);
            dNz(:,13)=0.25*a1.*a2.*d;
            dNz(:,14)=0.25*b1.*a2.*d;
            dNz(:,15)=0.25*b1.*b2.*d;
            dNz(:,16)=0.25*a1.*b2.*d;
            
        end
        function [node,dir]=SideDef(sidenum)
            nodes=[1,2,6,5,9,14,17,13
                2,3,7,6,10,15,18,14
                3,4,8,7,11,16,19,15
                1,5,8,4,13,20,16,12
                1,4,3,2,12,11,10,9
                5,6,7,8,17,18,19,20];
            dirs=[-2 1 2 -1 -3 3];
            node=nodes(sidenum,:);
            dir=dirs(sidenum)';
        end
        function [v,conn,elemtype]=GenSubCells(npts,type)
            if nargin<2,
                type='s';
            end
            tol=1e-14;
            [v,conn,elemtype]=ShapeFnc.Hex8.GenSubCells(npts,type);
            nv=size(v,1);
            v=[v
                (v(conn(:,2),:)+v(conn(:,1),:))/2
                (v(conn(:,3),:)+v(conn(:,2),:))/2
                (v(conn(:,4),:)+v(conn(:,3),:))/2
                (v(conn(:,4),:)+v(conn(:,1),:))/2
                (v(conn(:,5),:)+v(conn(:,1),:))/2
                (v(conn(:,6),:)+v(conn(:,2),:))/2
                (v(conn(:,7),:)+v(conn(:,3),:))/2
                (v(conn(:,8),:)+v(conn(:,4),:))/2
                (v(conn(:,5),:)+v(conn(:,6),:))/2
                (v(conn(:,6),:)+v(conn(:,7),:))/2
                (v(conn(:,7),:)+v(conn(:,8),:))/2
                (v(conn(:,8),:)+v(conn(:,5),:))/2];
            nc=size(conn,1);
            conn=[conn ...
                (nv+(1:nc))'...
                (nv+((nc+1):(2*nc)))' ...
                (nv+((2*nc+1):(3*nc)))' ...
                (nv+((3*nc+1):(4*nc)))' ...
                (nv+((4*nc+1):(5*nc)))' ...
                (nv+((5*nc+1):(6*nc)))' ...
                (nv+((6*nc+1):(7*nc)))' ...
                (nv+((7*nc+1):(8*nc)))' ...
                (nv+((8*nc+1):(9*nc)))' ...
                (nv+((9*nc+1):(10*nc)))' ...
                (nv+((10*nc+1):(11*nc)))' ...
                (nv+((11*nc+1):(12*nc)))'];
            [v,I,J]=unique(tol*round(v/tol),'rows');
            conn=J(conn);
            if npts==1,
                conn=conn';
            end
        end
        %%%%%
        function conn=ReverseElement(conn)
            conn=conn(:,[1 4 3 2 5 8 7 6 12 11 10 9 13 16 15 14 20 19 18 17]);
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
                        vi=1:20;
                    case 2
                        vi=[3 2 6 7 4 1 5 8 10 14 18 15 11 9 17 19 12 13 20 16];
                    case 3
                        vi=[1 5 6 2 4 8 7 3 13 17 14 9 12 20 18 10 16 19 15 11];
                end
                conn(i,:)=conn(i,vi);
                if vi2==2,
                    vi1=vi1(4:-1:1);
                    vi=[vi1+4 vi1 vi1([2 3 4 1])+16 vi1+12 vi1([2 3 4 1])+8];
                else
                    vi=[vi1 vi1+4 vi1+8 vi1+12 vi1+16];
                end
                conn(i,:)=conn(i,vi);
            end
        end
    end
end
