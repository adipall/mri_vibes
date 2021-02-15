classdef Shape2D < ShapeFnc.Shape3D
    properties
    end
    
    methods
        function obj=Shape2D(conn,nodexyz)
            obj=obj@ShapeFnc.Shape3D(conn,nodexyz);
            obj.Dimension=size(nodexyz,2);
        end
        function [pts]=Global2Local(obj,PTS,elemnum)
            switch obj.Dimension
                case 2
                    [pts]=Global2Local_2D(obj,PTS,elemnum);
                case 3
                    %error('Global2Local_3D: Not yet implemented')
                    [pts]=Global2Local_3D(obj,PTS,elemnum);
                otherwise 
                    error('Unreasonable dimension: dim = %d',obj.Dimension);
            end
        end
        function [pts]=Global2Local_2D(obj,PTS,elemnum)
            pts=obj.LocalCentroid;
            pts=ones(size(PTS,1),1)*pts(:,1:2);
            xiterp=1e16*ones(size(PTS,1),2);
            % force zeros in the z's
            err=1e16*ones(size(PTS,1),1);
            niter=1;
            Xc=obj.Coordinates(:,1);
            Yc=obj.Coordinates(:,2);
            
            %
            while (niter<=obj.MaxIter && any(err>obj.MaxError))
                N=obj.Interpolate(pts);
                if length(elemnum)>1,
                    f=PTS(:,1:2)-[sum(N.*Xc(obj.Connectivity(elemnum,:)),2) sum(N.*Yc(obj.Connectivity(elemnum,:)),2)];
                else
                    f=PTS(:,1:2)-[sum(N.*Xc(obj.Connectivity(elemnum,:))',2) sum(N.*Yc(obj.Connectivity(elemnum,:))',2)];
                end
                J1=obj.Jacobian(pts,elemnum);
                %% actually need the transpose of J1
                pts=pts+VMath.Solve22(J1(:,[1 3 2 4]),f); 
                
                niter=niter+1;
                err=sum((pts-xiterp).*(pts-xiterp),2);
                
                xiterp=pts;
            end
            if niter>obj.MaxIter,
                % if it can't converge in maxiter iterations, then set the
                % coordinates to nan and get on with life
                %
                %idx=find(err>obj.MaxError);
                %pts(idx,:)=repmat([NaN NaN NaN],length(idx),1);
                error('Maximum iterations exceeded')
            end
        end
        function [pts]=Global2Local_3D(obj,PTS,elemnum)
            pts=obj.LocalCentroid;
            pts=ones(size(PTS,1),1)*pts;
            xiterp=1e16*ones(size(PTS,1),2);
            % force zeros in the z's
            err=1e16*ones(size(PTS,1),1);
            niter=1;
            Xc=obj.Coordinates(:,1);
            Yc=obj.Coordinates(:,2);
            Zc=obj.Coordinates(:,3);
            %
            while (niter<=obj.MaxIter && any(err>obj.MaxError))
                N=obj.Interpolate(pts);
                if length(elemnum)>1,
                    f=PTS-[sum(N.*Xc(obj.Connectivity(elemnum,:)),2) sum(N.*Yc(obj.Connectivity(elemnum,:)),2) sum(N.*Zc(obj.Connectivity(elemnum,:)),2)];
                else
                    f=PTS-[sum(N.*Xc(obj.Connectivity(elemnum,:))',2) sum(N.*Yc(obj.Connectivity(elemnum,:))',2) sum(N.*Zc(obj.Connectivity(elemnum,:))',2)];
                end
                J1=obj.Jacobian(pts,elemnum);
                %% VMath stores non rectangular matrices kind of different
               
                pts=pts+VMath.PInv32(J1,f); 
                
                niter=niter+1;
                err=sum((pts-xiterp).*(pts-xiterp),2);
                
                xiterp=pts;
            end
            if niter>obj.MaxIter,
                % if it can't converge in maxiter iterations, then set the
                % coordinates to nan and get on with life
                %
                %idx=find(err>obj.MaxError);
                %pts(idx,:)=repmat([NaN NaN NaN],length(idx),1);
                error('Maximum iterations exceeded')
            end
        end
        %%%
        
        %%%
        function [invJ,detJ]=invJacobian(obj,pts,elemnum)
            J1=obj.Jacobian(pts,elemnum);
            [invJ,detJ]=VMath.Inv22(J1);
        end
        function [J1,dNx,dNy]=Jacobian(obj,pts,elemnum)
            dNx=obj.lderivX(pts);
            dNy=obj.lderivY(pts);

            x=reshape(obj.Coordinates(obj.Connectivity(elemnum,:),1),length(elemnum),size(obj.Connectivity,2));
            y=reshape(obj.Coordinates(obj.Connectivity(elemnum,:),2),length(elemnum),size(obj.Connectivity,2));
            
            switch obj.Dimension
                case 2
                    J1=[sum(x.*dNx,2) sum(y.*dNx,2) ...
                        sum(x.*dNy,2) sum(y.*dNy,2)];
                case 3
                    z=reshape(obj.Coordinates(obj.Connectivity(elemnum,:),3),length(elemnum),size(obj.Connectivity,2));
                    J1=[sum(x.*dNx,2) sum(y.*dNx,2) sum(z.*dNx,2) ...
                        sum(x.*dNy,2) sum(y.*dNy,2) sum(z.*dNy,2)];
                otherwise
                    error('Unknown dimension: %d',obj.Dimension)
            end
        end
        function [dNx,dNy,detJ]=Derivatives(obj,pts,elemnum)
            
            [J1,dNx,dNy]=obj.Jacobian(pts,elemnum);
           
            [dNx,dNy,detJ]=ShapeFnc.Shape2D.calcGlobalDerivatives(dNx,dNy,J1);
            
        end
        
        function [detJ]=calcDetJ(obj,pts,elemnum)
            J1=obj.Jacobian(pts,elemnum);
            detJ=J1(:,1).*J1(:,4)-J1(:,2).*J1(:,3);
        end
        
        function [J1,dAx,dAy]=SurfJacobian(obj,pts,sidenum,elemnum)
            if nargin==3,
                elemnum=1:size(obj.Connectivity,1);
            end
            
            [nodes,sdir]=obj.SideDef(sidenum);
            
            J1=obj.Jacobian(pts,elemnum);

            dAx=zeros(length(elemnum),1);
            dAy=zeros(length(elemnum),1);
            switch abs(sdir),
                case 1
                    dAx=J1(:,3);
                    dAy=J1(:,4);
                case 2
                    dAx=J1(:,1);
                    dAy=J1(:,2);
                case 3  %% this should only be used by triangles!!
                    dAx=J1(:,1)-J1(:,3);
                    dAy=J1(:,2)-J1(:,4);
                otherwise
                    error('Direction unknown %d in Shape2D.SurfJacobian',abs(sdir));
            end
            J1=sqrt(dAx.*dAx+dAy.*dAy);
            %           dAx=sign(sdir).*dAx;
            %           dAy=sign(sdir).*dAy;
        end
        function [B,a,Y]=QuadApprox(obj,elemnum)
            %% This generates a quadratic approximation to the volume of
            %% the element according to
            %% Rashid, International Journal of Numerical Methods in
            %% Engineering, vol 55, pp 431-450, 2002.
            if nargin==1,
                elemnum=1:size(obj.Connectivity,1);
            end
            %
            nintpts=obj.NumSurfIntPts+1;
            nelem=length(elemnum);
            %
            xx=zeros(nelem,3);
            a=zeros(nelem,2);
            A=zeros(nelem,1);
            x=obj.Coordinates(:,1);
            y=obj.Coordinates(:,2);
            elem=repmat(elemnum(:)',nintpts,1);
            elem=elem(:);
            for i=1:obj.NumSides,
                [pts,we]=obj.getSurfIntPts(i,nintpts);
                sides=obj.SideDef(i);
                
                sidnod=obj.Connectivity(elemnum,sides);
                
                N=obj.Interpolate(pts);
                J1=obj.SurfJacobian(repmat(pts,nelem,1),i,elem);
                J1r=reshape(repmat(we,nelem,1).*J1,nintpts,nelem)';
                %%%%%%
                
                a(:,1)=a(:,1)+sum(J1r.*(x(sidnod)*N(:,sides)'),2);
                a(:,2)=a(:,2)+sum(J1r.*(y(sidnod)*N(:,sides)'),2);
                A=A+sum(J1r,2);
            end
            a=a./A(:,[1 1]);
            ns=length(sides);
            for i=1:obj.NumSides,
                [pts,we]=obj.getSurfIntPts(i,nintpts);
                sides=obj.SideDef(i);
                sidnod=obj.Connectivity(elemnum,sides);
                N=obj.Interpolate(pts);
                J1=obj.SurfJacobian(repmat(pts,nelem,1),i,elem);
                J1r=reshape(repmat(we,nelem,1).*J1,nintpts,nelem)';
                %%%%%%
                xx(:,1)=xx(:,1)+sum(J1r.*(((x(sidnod)-a(:,ones(ns,1)))*N(:,sides)').*((x(sidnod)-a(:,ones(ns,1)))*N(:,sides)')),2);
                xx(:,3)=xx(:,3)+sum(J1r.*(((x(sidnod)-a(:,ones(ns,1)))*N(:,sides)').*((y(sidnod)-a(:,2*ones(ns,1)))*N(:,sides)')),2);
                xx(:,2)=xx(:,2)+sum(J1r.*(((y(sidnod)-a(:,2*ones(ns,1)))*N(:,sides)').*((y(sidnod)-a(:,2*ones(ns,1)))*N(:,sides)')),2);
            end
            A=A./(xx(:,1).*xx(:,2)-xx(:,3).*xx(:,3));
            B=[xx(:,2) xx(:,1) -xx(:,3)].*A(:,[1 1 1]);
            if nargout>2,
                nintpts=obj.NumIntPts*obj.NumIntPts;
                elem=repmat(elemnum(:)',nintpts,1);
                elem=elem(:);

                [pts,we]=obj.getIntPts;
                pt=repmat(pts,nelem,1);
                N=obj.Interpolate(pts);

                J1=reshape(obj.calcDetJ(pt,elem),nintpts,nelem)';
                Y=zeros(nelem,2);
                Y(:,1)=((x(obj.Connectivity(elemnum,:))*N').*J1)*we;
                Y(:,2)=((y(obj.Connectivity(elemnum,:))*N').*J1)*we;
                
            end
        end
        %%%%
    end
    
    methods (Static)
        function [dNx,dNy,detJ]=calcGlobalDerivatives(dNx,dNy,J1)
            
            f=zeros(size(dNx,1),2,size(dNx,2));
            f(:,1,:)=dNx;
            f(:,2,:)=dNy;
            
            [sol,detJ]=VMath.Solve22(J1,f);
            dNx=squeeze(sol(:,1,:));
            dNy=squeeze(sol(:,2,:));
                
        end
    end
end
