classdef Tet10 < ShapeFnc.Shape3D
   properties
       PltConn = [1 5 2 6 3 7 1 4 10 3 6 2 9 4]
   end

   methods
      %% All Interp matrices are nPoints x nNodes
      function obj=Tet10(conn,nodexyz)
         obj=obj@ShapeFnc.Shape3D(conn,nodexyz);
         obj.NumSides=4;
         obj.NumIntPts=2;
         obj.NumMassIntPts=obj.NumIntPts;
         obj.NumSurfIntPts=2;
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
              [pts,weights]=ShapeFnc.Shape3D.TetIntegration(n);
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
   end
   methods (Static)
      function N=Interpolate(pts)
         % corner nodes  the z direction needs to go first
         z=1-pts(:,1)-pts(:,2)-pts(:,3);
         N=[(2*z-1).*z ...
             (2*pts(:,1)-1).*pts(:,1) ... 
             (2*pts(:,2)-1).*pts(:,2) ...
             (2*pts(:,3)-1).*pts(:,3)  ...
             4*pts(:,1).*z ...         %8
             4*pts(:,1).*pts(:,2) ...  %5
             4*pts(:,2).*z ...         %9
             4*pts(:,3).*z ...         % 10
             4*pts(:,1).*pts(:,3) ...  %7
             4*pts(:,2).*pts(:,3)];    %6
      end
      function [pts]=LocalCentroid
          pts=[1/4 1/4 1/4];
      end
      function [xb,yb,zb]=LocalBounds
            xb=[0 1];
            yb=[0 1];
            zb=[0 1];
      end
      %%
  
      function dNx=lderivX(pts)
          n=size(pts,1);
          z=1-pts(:,1)-pts(:,2)-pts(:,3);
          dNx=[(1-4*z) 4*pts(:,1)-1 zeros(n,1) zeros(n,1) ...
              4*(z-pts(:,1)) 4*pts(:,2) -4*pts(:,2)...
               -4*pts(:,3) 4*pts(:,3) zeros(n,1)];
      end
      function dNy=lderivY(pts)
          n=size(pts,1);
          z=1-pts(:,1)-pts(:,2)-pts(:,3);
          dNy=[1-4*z zeros(n,1) 4*pts(:,2)-1 zeros(n,1) ...
              -4*pts(:,1) 4*pts(:,1) 4*(z-pts(:,2)) ...
              -4*pts(:,3) zeros(n,1) 4*pts(:,3)];
      end
      function dNz=lderivZ(pts)
          n=size(pts,1);
          z=1-pts(:,1)-pts(:,2)-pts(:,3);
          dNz=[1-4*z zeros(n,1) zeros(n,1) 4*pts(:,3)-1 ...
               -4*pts(:,1) zeros(n,1) -4*pts(:,2) 4*(z-pts(:,3)) ...
               4*pts(:,1) 4*pts(:,2)];
      end
      function [node,dir]=SideDef(sidenum)
          nodes=[1,2,4,5,9,8
              2,3,4,6,10,9
              1,4,3,8,10,7
              1,3,2,7,6,5];
          dirs=[3 1 2 4];
          node=nodes(sidenum,:);
          dir=dirs(sidenum)';
      end
      function [v,conn,elemtype]=GenSubCells(npts,type)
          if type==1,
              v=[];
              conn=1:10;
              elemtype='tet';
          else
              
              error('double mesh not implemented for tets')
          end
      end
      %%%%%
      function conn=ReverseElement(conn)
            conn=conn(:,[1 3 2 4 7 6 5 8 10 9]);
      end
      function conn=RandomizeConn(conn)
          m1=[1 2 3 4 5 6 7 8 9 10
              2 3 1 4 6 7 5 9 10 8
              3 1 2 4 7 5 6 10 8 9
              %%
              2 4 3 1 9 10 6 5 8 7
              3 2 4 1 6 9 10 7 5 8
              4 3 2 1 10 6 9 8 7 5
              %%
              1 3 4 2 7 10 8 5 6 9
              4 1 3 2 8 7 10 9 5 6
              3 4 1 2 10 8 7 6 9 5];
          for i=1:size(conn,1),
              vi1=m1(randi(9,1),:);
              
              conn(i,:)=conn(i,vi1);
          end
      end
   end
end
