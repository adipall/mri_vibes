classdef Wed15 < ShapeFnc.Shape3D
   properties
       PltConn = [1 7 2 8 3 9 1 10 4 13 5 14 6 15 4 15 6 12 3 8 2 11 5]
   end

   methods
      %% All Interp matrices are nPoints x nNodes
      function obj=Wed15(conn,nodexyz)
         obj=obj@ShapeFnc.Shape3D(conn,nodexyz);
         obj.NumSides=5;
         obj.NumIntPts=4;
         obj.NumMassIntPts=obj.NumIntPts+1;
         obj.NumSurfIntPts=4;
      end
      function [pts,weights,idnum]=getIntPts(obj,n)
          idnum=[];
          if nargin<2,
              n=obj.NumIntPts;
          end
          
          [pts,weights,idnum]=ShapeFnc.Shape3D.WedIntegration(n);
      end
      function [pts,weights,idnum]=getIntPts_old(obj,n)
          idnum=[];
          if nargin<2,
              n=obj.NumIntPts;
          end
          
          [pt,wt]=ShapeFnc.Shape3D.TriIntegration(n);
          [pz,wz]=ShapeFnc.Shape3D.GaussWeights(n);
          pz=repmat(pz',size(pt,1),1);
          wz=repmat(wz',size(pt,1),1);
          pts=[repmat(pt,n,1) pz(:)];
          weights=wz(:).*repmat(wt,n,1);
          if nargout==3,
              idnum=ShapeFnc.Shape3D.genGaussPointNames(pts);
          end
          
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
                      case 2
                          pts=[pt(:) 1-pt(:) pz(:)];
                      case 3
                          pts=[zeros(numel(pt),1) pt(:) pz(:)];
                      case 1
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
                  error('Unknown side in Wed15')
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
   end
   methods (Static)
      function N=Interpolate(pts)
         % corner nodes
         N=zeros(size(pts,1),15);
         a1=(1-pts(:,3));
         b1=(1+pts(:,3));
         d=(1-pts(:,3).*pts(:,3));
         t=1-pts(:,1)-pts(:,2);
         N(:,1:6)=0.5*[t.*(2*t-1).*a1 - t.*d ...
             pts(:,1).*(2*pts(:,1)-1).*a1 - pts(:,1).*d ...
             pts(:,2).*(2*pts(:,2)-1).*a1 - pts(:,2).*d ...
             t.*(2*t-1).*b1 - t.*d ...
             pts(:,1).*(2*pts(:,1)-1).*b1 - pts(:,1).*d ...
             pts(:,2).*(2*pts(:,2)-1).*b1 - pts(:,2).*d];
         %top ones midsides of triangles
         N(:,7)=2*pts(:,1).*t.*a1;
         N(:,8)=2*pts(:,1).*pts(:,2).*a1;
         N(:,9)=2*pts(:,2).*t.*a1;
         % mi
         N(:,10)=t.*d;
         N(:,11)=pts(:,1).*d;
         N(:,12)=pts(:,2).*d;
         % bottom ones midsides of triangles
         N(:,13)=2*pts(:,1).*t.*b1;
         N(:,14)=2*pts(:,1).*pts(:,2).*b1;
         N(:,15)=2*pts(:,2).*t.*b1;
      end
      function [pts]=LocalCentroid
         pts=[1/3 1/3 0];
      end
      function pts=LocalNode
          pts=[0  0   -1
              1   0   -1
              0   1   -1
              0   0    1
              1   0    1
              0   1    1
              0.5 0   -1
              0.5 0.5 -1
              0   0.5 -1
              0   0    0
              1   0    0
              0   1    0
              0.5 0    1
              0.5 0.5  1
              0   0.5  1];
      end
      function [xb,yb,zb]=LocalBounds
            xb=[0 1];
            yb=[0 1];
            zb=[-1 1];
      end

      function dNx=lderivX(pts)
          dNx=zeros(size(pts,1),15);
          a1=(1-pts(:,3));
          b1=(1+pts(:,3));
          d=(1-pts(:,3).*pts(:,3));
          t=1-pts(:,1)-pts(:,2);
          z=zeros(size(pts,1),1);
          dNx(:,1:6)=0.5*[(1-4*t).*a1 + d ... %dN/dt*dt/dr
              (4*pts(:,1) -1).*a1 - d ...
              z ...
              (1-4*t).*b1 + d ...
              (4*pts(:,1) -1).*b1 - d ...
              z];
          %
          dNx(:,7:15)=[2*(t-pts(:,1)).*a1 ...
              2*pts(:,2).*a1 ...
              -2*pts(:,2).*a1 ...
              -d ...
              d ...
              z ...
              2*(t-pts(:,1)).*b1 ...
              2*pts(:,2).*b1 ...
              -2*pts(:,2).*b1];
      end
      function dNy=lderivY(pts)
          dNy=zeros(size(pts,1),15);
          a1=(1-pts(:,3));
          b1=(1+pts(:,3));
          d=(1-pts(:,3).*pts(:,3));
          t=1-pts(:,1)-pts(:,2);
          z=zeros(size(pts,1),1);
          dNy(:,1:6)=0.5*[(1-4*t).*a1 + d  ...
              z ...
              (4*pts(:,2) -1).*a1 - d ...
              (1-4*t).*b1 + d ...
              z ...
              (4*pts(:,2) -1).*b1 - d];
          dNy(:,7:15)=[-2*pts(:,1).*a1 ...
              2*pts(:,1).*a1 ...
              2*(t-pts(:,2)).*a1 ...
              -d ...
              z ...
              d ...
              -2*pts(:,1).*b1 ...
              2*pts(:,1).*b1 ...
              2*(t-pts(:,2)).*b1];
      end
      function dNz=lderivZ(pts)
          d=-2*pts(:,3);
          t=1-pts(:,1)-pts(:,2);
          dNz=zeros(size(pts,1),15);
          dNz(:,1:6)=0.5*[-t.*(2*t-1) - t.*d ...
              -pts(:,1).*(2*pts(:,1)-1) - pts(:,1).*d ...
              -pts(:,2).*(2*pts(:,2)-1) - pts(:,2).*d ...
              t.*(2*t-1) - t.*d ...
              pts(:,1).*(2*pts(:,1)-1) - pts(:,1).*d ...
              pts(:,2).*(2*pts(:,2)-1) - pts(:,2).*d];
          dNz(:,7:15)=[-2*pts(:,1).*t ...
              -2*pts(:,1).*pts(:,2) ...
              -2*pts(:,2).*t ...
              t.*d ...
              pts(:,1).*d ...
              pts(:,2).*d ...
              2*pts(:,1).*t ...
              2*pts(:,1).*pts(:,2) ...
              2*pts(:,2).*t];
      end
      function [node,dir]=SideDef(sidenum)
          switch sidenum
              case 1
                  node=[1,2,5,4,7,11,13,10];
              case 2
                  node=[2,3,6,5,8,12,14,11];
              case 3
                  node=[1,4,6,3,10,15,12,9];
              case 4
                  node=[1,3,2,9,8,7];
              case 5
                  node=[4,5,6,13,14,15]; % 3 node side
          end
          dirs=[2 4 1 -3 -3];
          dir=dirs(sidenum)';
      end
      %%%%%%
      function [v,conn,elemtype]=GenSubCells(npts,type)
            if nargin<2,
                type='s';
            end
            tol=1e-14;
            [v,conn,elemtype]=ShapeFnc.Wed6.GenSubCells(npts,type);
            nv=size(v,1);
            v=[v
                (v(conn(:,2),:)+v(conn(:,1),:))/2
                (v(conn(:,3),:)+v(conn(:,2),:))/2
                (v(conn(:,3),:)+v(conn(:,1),:))/2
                (v(conn(:,4),:)+v(conn(:,1),:))/2
                (v(conn(:,5),:)+v(conn(:,2),:))/2
                (v(conn(:,6),:)+v(conn(:,3),:))/2
                (v(conn(:,4),:)+v(conn(:,5),:))/2
                (v(conn(:,5),:)+v(conn(:,6),:))/2
                (v(conn(:,6),:)+v(conn(:,4),:))/2];
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
                (nv+((8*nc+1):(9*nc)))'];
            [v,I,J]=unique(tol*round(v/tol),'rows');
            conn=J(conn);
            if npts==1,
                conn=conn';
            end
      end
      %%%%%%%
      function conn=ReverseElement(conn)
            conn=conn(:,[1 3 2 4 6 5 9 8 7 10 12 11 15 14 13]);
      end
      function conn=RandomizeConn(conn)
            m1=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
                2 3 1 5 6 4 8 9 7 11 12 10 14 15 13
                3 1 2 6 4 5 9 7 8 12 10 11 15 13 14
                4 6 5 1 3 2 15 14 13 10 12 11 9 8 7
                6 5 4 3 2 1 14 13 15 12 11 10 8 7 9
                5 4 6 2 1 3 13 15 14 11 10 12 7 9 8];
            for i=1:size(conn,1),
                vi1=m1(randi(6,1),:);
                
                conn(i,:)=conn(i,vi1);
            end
        end
   end
end
