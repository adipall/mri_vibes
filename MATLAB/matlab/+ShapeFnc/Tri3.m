classdef Tri3 < ShapeFnc.Shape2D
   properties 
       PltConn = [1 2 3 1]
   end

   methods
      %% All Interp matrices are nPoints x nNodes
      function obj=Tri3(conn,nodexyz)
         obj=obj@ShapeFnc.Shape2D(conn,nodexyz);
         obj.NumSides=3;
         obj.NumIntPts=2;
         obj.NumMassIntPts=obj.NumIntPts;
         obj.NumSurfIntPts=2;
      end
      function [pts]=Global2Local2(obj,nodexyz,elemnum)
          %% this isn't robust
          %
          %%x*X1+x*X2+t*X3=X
          %%x*Y1+x*Y2+t*Y3=Y
          %%x+y+t=1 => t=1-x-y;
          
          Xc=obj.Coordinates(:,1);
          Yc=obj.Coordinates(:,2);
          Xc=Xc(obj.Connectivity(elemnum,:));
          Yc=Yc(obj.Connectivity(elemnum,:));
          %%
          a1=(nodexyz(:,1)-Xc(:,1))./(Xc(:,1)-Xc(:,3));
          a2=(Xc(:,2)-Xc(:,3))./(Xc(:,1)-Xc(:,3));
          %%
          b1=(nodexyz(:,2)-Yc(:,3))./(Yc(:,1)-Yc(:,3));
          b2=(Yc(:,2)-Yc(:,3))./(Yc(:,1)-Yc(:,3));
          
          pts=zeros(size(nodexyz,1),3);
          pts(:,2)=(b1-a1)./(b2-a2);
          pts(:,1)=a1-pts(:,2).*a2;
      end
      function [vol]=Volume(obj,elements)
          if nargin==1,
              elements=1:size(obj.Connectivity,1);
          end
          
          x=obj.Coordinates(:,1);
          y=obj.Coordinates(:,2);
          conn=obj.Connectivity(elements,:);
          vol=0.5*((x(conn(:,2)).*y(conn(:,3))-y(conn(:,2)).*x(conn(:,3))) ...
                        -x(conn(:,1)).*(y(conn(:,3))-y(conn(:,2))) ...
                        +y(conn(:,1)).*(x(conn(:,3))-x(conn(:,2))));
      end
      function [pts,weights,idnum]=getIntPts(obj,n)
          if nargin==1,
              n=obj.NumIntPts;
          end
          if length(n)==1,
              [pts,weights]=ShapeFnc.Shape3D.TriIntegration(n);
          else
              error('Incorrect number of gauss points specified')
          end
          if nargout==3,
              idnum=ShapeFnc.Shape3D.genGaussPointNames(pts);
          end
      end
      function [pts,weights]=getGenIntPts(obj,n,limits)
          if nargin==1,
              n=obj.NumIntPts;
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
                (c(conn(:,1),:)+c(conn(:,3),:))]/2;
            n=[c(conn(:,2),:)-c(conn(:,1),:) ...
                c(conn(:,3),:)-c(conn(:,2),:) ...
                c(conn(:,1),:)-c(conn(:,3),:)];
            nn=[sqrt(sum(n(:,[1 2]).*n(:,[1 2]),2)) ...
                sqrt(sum(n(:,[3 4]).*n(:,[3 4]),2)) ...
                sqrt(sum(n(:,[5 6]).*n(:,[5 6]),2))];
            n=[n(:,2)./nn(:,1) -n(:,1)./nn(:,1) ...
                n(:,4)./nn(:,2) -n(:,3)./nn(:,2) ...
                n(:,6)./nn(:,3) -n(:,5)./nn(:,3)];

            r=rc-[PTS PTS PTS];
            r=r.*n;
            num=[r(:,1)+r(:,2) r(:,3)+r(:,4) r(:,5)+r(:,6)];
            t=n.*[t t t];
            den=[t(:,1)+t(:,2) t(:,3)+t(:,4) t(:,5)+t(:,6)];
            
            alpha=1./sort(den./num,2);
            switch 1
                case 1
                    %%for particle tracing
                    vi=3;
                case 2
                    %% for contact
                    vi=1;
            end
                    
            PTSI=PTS+[alpha(:,vi) alpha(:,vi)].*vel./mag(:,[1 1]);
           
        end
   end
   methods (Static)
      function N=Interpolate(pts)
          % use area coordinates
          % corner nodes
          l3=1-pts(:,1)-pts(:,2);
          N=[pts(:,1) pts(:,2) l3];
      end
      function [pts]=LocalCentroid
         pts=[1/3 1/3];
      end
      function [xb,yb,zb]=LocalBounds
            xb=[0 1];
            yb=[0 1];
            zb=[0 0];
      end
      %%%
      function dNx=lderivX(pts)
          n=size(pts,1);
          dNx=[ones(n,1) zeros(n,1) -ones(n,1)];
      end
      function dNy=lderivY(pts)
          n=size(pts,1);
          dNy=[zeros(n,1) ones(n,1) -ones(n,1) ];
      end
      function dNz=lderivZ(pts)
          dNz=zeros(size(pts,1),3);
      end
      %%%
      function [node,dir]=SideDef(sidenum)
          nodes=[1,2
              2,3
              3,1];
          dirs=[3 1 2];
          node=nodes(sidenum,:);
          dir=dirs(sidenum)';     
      end
      function [v,conn,elemtype]=GenSubCells(npts,type)
          if nargin<2,
              type='s';
          end
          if strcmpi(type,'g')
              if npts>5,
                  error('Tri.GenSubCells not defined for gauss points greater than 5th order')
              end
          end
          [v,conn,elemtype]=ShapeFnc.Tri3.dbl_mesh(npts,type);
      end
      function [v,conn,elemtype]=dbl_mesh(npts,type)
          if strcmpi(type,'s')
              elemtype='tri';
              switch npts
                  case 1
                      v=[1 0;0 1;0 0];
                      conn=[1 2 3];
                  case 2
                      v=[1 0;0 1;0 0;.5 .5;0 .5;.5 0];
                      conn=[1 4 6;
                          6 4 5;
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
          elseif strcmpi(type,'g'),
              switch npts
                  case 1
                      v=[1 0;0 1;0 0];
                      conn=[1 2 3];
                      elemtype='tri';
                  case 2
                      v=[1 0;0 1;0 0;.5 .5;0 .5];
                      conn=[1 5 3;1 4 5;4 2 5];
                      elemtype='tri';
                  case 3
                      v=[];
                      conn=[];
                      elemtype=[];
                  case 4
                      v=[1 0;0 1;0 0;1/3 1/3];
                      conn=[1 2 4;3 4 2;1 4 3];
                      elemtype='tri';
                  case 5
                      v=[1 0;0 1;0 0;1/3 0;2/3 0;0 1/3;0 2/3;.5 .5];
                      conn=[1 8 5;7 8 2;3 4 6;5 7 6;5 6 4;5 8 7];
                      elemtype='tri';
                  otherwise
                      error('Unknown number of points')
              end
          else
              error('Unknown type')
          end
      end
      function conn=ReverseElement(conn)
          conn=conn(:,[1 3 2]);
      end
      function conn=RandomizeConn(conn)
          m1=[1 2 3
              2 3 1
              3 1 2];
          for i=1:size(conn,1),
              vi1=m1(randi(3,1),:);
              
              conn(i,:)=conn(i,vi1);
          end
      end
   end
end
