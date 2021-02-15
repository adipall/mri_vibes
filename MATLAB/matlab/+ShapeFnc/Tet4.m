classdef Tet4 < ShapeFnc.Shape3D
   properties
       PltConn = [1 2 3 1 4 3 2 4]
   end

   methods
      %% All Interp matrices are nPoints x nNodes
      function obj=Tet4(conn,nodexyz)
         obj=obj@ShapeFnc.Shape3D(conn,nodexyz);
         obj.NumSides=4;
         obj.NumIntPts=1;
         obj.NumMassIntPts=obj.NumIntPts;
         obj.NumSurfIntPts=1;
      end
      function [pts,weights,idnum]=getIntPts(obj,n)
          if nargin==1,
              n=obj.NumIntPts;
          end
          if length(n)==1,
              [pts,weights,idnum]=ShapeFnc.Shape3D.TetIntegration(n);
          else
              error('Incorrect number of gauss points specified')
          end 
      end
      function [pts,weights]=getGenIntPts(obj,n,limits)
          if nargin==1,
              n=obj.NumIntPts;
          end
          if length(n)==1,
              [pts,weights]=Shape3D.TetIntegration(n);
          elseif length(n)==3,
              [pts,weights]=Shape3D.TetIntegration(n);
          else
              error('Incorrect number of gauss points specified')
          end
      end
      
      function [pts,weights]=getSurfIntPts(obj,side,n)
          [node,dir]=obj.SideDef(side);
          if nargin==2,
              n=obj.NumIntPts;
          end
          if length(n)==1,
              [pts,weights]=ShapeFnc.Shape3D.TriIntegration(n);
          else
              error('Incorrect number of gauss points specified')
          end
          switch side
              case 1 % correct
                  pts=[pts(:,1) zeros(size(pts,1),1) pts(:,2)];
              case 2
                  pts=[pts 1-sum(pts,2)];
              case 3 % correct
                  pts=[zeros(size(pts,1),1) pts];
              case 4 % correct
                  pts=[pts zeros(size(pts,1),1)];
              otherwise
                  error('Unknown direction in Shape3D.')
          end
      end
      function [pts]=Global2Local(obj,nodexyz,elemnum)
          
          %%x*X1+y*X2+z*X3+t*X4=X
          %%x*Y1+y*Y2+z*Y3+t*Y4=Y
          %%x*Z1+y*Z2+z*Y3+t*Y4=Z
          %%x+y+z+t=1 => t=1-x-y-z;
          Xc=obj.Coordinates(:,1);
          Yc=obj.Coordinates(:,2);
          Zc=obj.Coordinates(:,3);
          
          Xc=Xc(obj.Connectivity(elemnum,:));
          Yc=Yc(obj.Connectivity(elemnum,:));
          Zc=Zc(obj.Connectivity(elemnum,:));
          %%
          %
          
          pts=ShapeFnc.Shape3D.inv33([(Xc(:,1)-Xc(:,4)) (Xc(:,2)-Xc(:,4)) (Xc(:,3)-Xc(:,4)) ...
              (Yc(:,1)-Yc(:,4)) (Yc(:,2)-Yc(:,4)) (Yc(:,3)-Yc(:,4)) ...
              (Zc(:,1)-Zc(:,4)) (Zc(:,2)-Zc(:,4)) (Zc(:,3)-Zc(:,4))], ...
              [(nodexyz(:,1)-Xc(:,4)) (nodexyz(:,2)-Yc(:,4)) (nodexyz(:,3)-Zc(:,4))]);
         
      end
      function [vol]=Volume(obj,elems)
          if nargin==1,
              elems=1:size(obj.Connectivity,1);
          end
          
          x=obj.Coordinates(:,1);
          y=obj.Coordinates(:,2);
          z=obj.Coordinates(:,3);
          
          L0=[x(obj.Connectivity(elems,2))-x(obj.Connectivity(elems,1)) ...
              y(obj.Connectivity(elems,2))-y(obj.Connectivity(elems,1))...
              z(obj.Connectivity(elems,2))-z(obj.Connectivity(elems,1))];
          %
          L2=[x(obj.Connectivity(elems,1))-x(obj.Connectivity(elems,3)) ...
              y(obj.Connectivity(elems,1))-y(obj.Connectivity(elems,3))...
              z(obj.Connectivity(elems,1))-z(obj.Connectivity(elems,3))];
          %
          L3=[x(obj.Connectivity(elems,4))-x(obj.Connectivity(elems,1)) ...
              y(obj.Connectivity(elems,4))-y(obj.Connectivity(elems,1))...
              z(obj.Connectivity(elems,4))-z(obj.Connectivity(elems,1))];
          
          vol=(L2(:,1).*(L0(:,2).*L3(:,3)-L0(:,3).*L3(:,2)) ...
              - L2(:,2).*(L0(:,1).*L3(:,3)-L0(:,3).*L3(:,1)) ...
              +L2(:,3).*(L0(:,1).*L3(:,2)-L0(:,2).*L3(:,1)))/6;
      end
      function [PTSI]=GlobalIntercept(obj,PTS,vel,elemnum)
          if nargin==3,
              elemnum=1:size(obj.Connectivity,1);
          end
          mag=sum(vel.*vel,2);
          t=vel./mag(:,[1 1]);
          c=obj.Coordinates(:,1:2);
          conn=obj.Connectivity(elemnum,:);
          rc=[(c(conn(:,1),:)+c(conn(:,2),:))+c(conn(:,4),:) ...
              (c(conn(:,2),:)+c(conn(:,3),:))+c(conn(:,4),:) ...
              (c(conn(:,1),:)+c(conn(:,4),:))+c(conn(:,3),:) ...
              (c(conn(:,1),:)+c(conn(:,3),:))+c(conn(:,2),:)]/2;
          n=[VMath.Cross(c(conn(:,4),:)-c(conn(:,1),:),c(conn(:,4),:)-c(conn(:,2),:)) ... %%% n will be the cross product of two edges on each face
              VMath.Cross(c(conn(:,4),:)-c(conn(:,2),:),c(conn(:,4),:)-c(conn(:,3),:)) ...
              VMath.Cross(c(conn(:,4),:)-c(conn(:,3),:),c(conn(:,4),:)-c(conn(:,1),:)) ...
              VMath.Cross(c(conn(:,3),:)-c(conn(:,2),:),c(conn(:,3),:)-c(conn(:,1),:))];
          nn=[sqrt(sum(n(:,[1 2 3]).*n(:,[1 2 3]),2)) ...
              sqrt(sum(n(:,[4 5 6]).*n(:,[4 5 6]),2)) ...
              sqrt(sum(n(:,[7 8 9]).*n(:,[7 8 9]),2)) ...
              sqrt(sum(n(:,[10 11 12]).*n(:,[10 11 12]),2))];
          n=[n(:,1:3)./nn(:,[1 1 1]) n(:,4:6)./nn(:,[2 2 2]) ...
              n(:,7:9)./nn(:,[3 3 3]) n(:,10:12)./nn(:,[4 4 4])];
          
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
      %%%
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
              (coord(conn(:,2),:)+coord(conn(:,4),:))/2
              (coord(conn(:,3),:)+coord(conn(:,4),:))/2];
          conn=[conn ((n+1):(n+nelem))' ((n+nelem+1):(n+2*nelem))' ...
              ((n+2*nelem+1):(n+3*nelem))' ((n+3*nelem+1):(n+4*nelem))' ...
              ((n+4*nelem+1):(n+5*nelem))' ((n+5*nelem+1):(n+6*nelem))'];
      end
   end
   methods (Static)
      function N=Interpolate(pts)
         % corner nodes evidently z has to be defined first
         z=1-pts(:,1)-pts(:,2)-pts(:,3);
         N=[z pts(:,1) pts(:,2) pts(:,3)];
      end
      
      function [pts]=LocalCentroid
          pts=[1/4 1/4 1/4];
      end
      function [xb,yb,zb]=LocalBounds
            xb=[0 1];
            yb=[0 1];
            zb=[0 1];
      end
      %%%

      function dNx=lderivX(pts)
          n=size(pts,1);
          dNx=[-ones(n,1) ones(n,1) zeros(n,1) zeros(n,1)];
      end
      function dNy=lderivY(pts)
          n=size(pts,1);
          dNy=[-ones(n,1) zeros(n,1)  ones(n,1) zeros(n,1)];
      end
      function dNz=lderivZ(pts)
          n=size(pts,1);
          dNz=[-ones(n,1) zeros(n,1) zeros(n,1) ones(n,1)];
      end
      function [node,dir]=SideDef(sidenum)
          nodes=[1,2,4
              2,3,4
              1,4,3
              1,3,2];
          dirs=[3 1 2 4];
          node=nodes(sidenum,:);
          dir=dirs(sidenum)';
      end
      %%%
      function [v,conn,elemtype]=GenSubCells(npts,type)
          if strcmpi(type,'g')
              if npts>2,
                  error('Tri.GenSubCells not defined for gauss points greater than 2 order')
              else
                v=[1 0 0;0 1 0;0 0 1;0 0 0];
                conn=[1 2 3 4];
                elemtype='tet';
              end
          else
              [v,conn]=ShapeFnc.Tet4.dbl_mesh(npts,type);
              elemtype='tet';
          end
          
      end
      function [v,conn]=dbl_mesh(npts,type)
          % base
          % 3
          % |\
          % 7-6
          % |  \
          % 1-5-2
          %
          % top
          %  9
          %  |\
          %  7-8
          if strcmpi(type,'g'),
              error('Gauss double mesh not implemented for Tets')
          end
          switch npts
              case 1
              v=[1 0 0;0 1 0;0 0 1;0 0 0];
              conn=[1 2 3 4];
              case 2
                  v=[1 0;0 1;0 0;.5 .5;0 .5;.5 0];
                  conn=[6 4 5;
                      1 4 6;
                      4 2 5;
                      6 5 3];
              case 3
                  v=[1 0;0 1;0 0;1/3 2/3;2/3 1/3;0 1/3;0 2/3;1/3 0;2/3 0;1/3 1/3];
                  conn=[3 8 6
                      6 10 7
                      7 4 2
                      4 7 10
                      4 10 5
                      8 9 10
                      6 8 10
                      5 10 9
                      1 5 9];
              case 4
                  v=[1 0;0 1;0 0;.25 .75;.5 .5;.75 .25;0 .75;0 .5;0 .25;.75 0;.5 0;.25 0; ...
                      .25 .25;.25 .5;.5 .25];
                  conn=[4 2 7
                      7 8 14
                      7 14 4
                      4 14 5
                      8 9 13
                      8 13 14
                      14 13 15
                      14 15 5
                      1 6 10
                      11 10 15
                      10 6 15
                      12 11 13
                      13 11 15
                      3 12 9
                      9 12 13
                      5 15 6];
              otherwise
                  error('Unknown number of points')
                      
          end
      end
      %%%%%
      function conn=ReverseElement(conn)
            conn=conn(:,[1 3 2 4]);
      end
      function conn=RandomizeConn(conn)
          m1=[1 2 3 4
              2 3 1 4
              3 1 2 4
              %%
              2 4 3 1
              3 2 4 1
              4 3 2 1
              %%
              1 3 4 2
              4 1 3 2
              3 4 1 2];
          for i=1:size(conn,1),
              vi1=m1(randi(9,1),:);
              
              conn(i,:)=conn(i,vi1);
          end
      end
   end
end
