classdef Shape1D < ShapeFnc.Shape3D
    properties
    end
    
    methods
        function obj=Shape1D(conn,nodexyz)
            obj=obj@ShapeFnc.Shape3D(conn,nodexyz);
            obj.Dimension=1;
        end
        function [pts]=Global2Local(obj,PTS,elemnum)
            pts=obj.LocalCentroid;
            pts=ones(size(PTS,1),1)*pts(:,1);
            xiterp=1e16*ones(size(PTS,1),1);
            % force zeros in the z's
            err=1e16*ones(size(PTS,1),1);
            niter=1;
            Xc=obj.Coordinates(:,1);
            
            %
            while (niter<=obj.MaxIter && any(err>obj.MaxError))
                N=obj.Interpolate(pts);
                if length(elemnum)>1,
                    f=PTS(:,1)-sum(N.*Xc(obj.Connectivity(elemnum,:)),2);
                else
                    f=PTS(:,1)-sum(N.*Xc(obj.Connectivity(elemnum,:))',2);
                end
                J1=obj.Jacobian(pts,elemnum);
                %% actually need the transpose of J1
                pts=pts+f./J1(:,1); 
                
                niter=niter+1;
                err=(pts-xiterp).*(pts-xiterp);
                
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
            detJ=obj.Jacobian(pts,elemnum);
            invJ=1./detJ;
        end
        function [J1,dNx]=Jacobian(obj,pts,elemnum)
            dNx=obj.lderivX(pts);
           
            x=reshape(obj.Coordinates(obj.Connectivity(elemnum,:),1),length(elemnum),size(obj.Connectivity,2));
            
            J1=sum(x.*dNx,2);
        end
        function [dNx,detJ]=Derivatives(obj,pts,elemnum)
            
            [J1,dNx]=obj.Jacobian(pts,elemnum);
           
            [dNx,detJ]=ShapeFnc.Shape1D.calcGlobalDerivatives(dNx,J1);
            
        end
        
        function [detJ]=calcDetJ(obj,pts,elemnum)
            J1=obj.Jacobian(pts,elemnum);
            detJ=J1(:,1);
        end
        
        function [J1,dAx]=SurfJacobian(obj,pts,sidenum,elemnum)
            if nargin==3,
                elemnum=1:size(obj.Connectivity,1);
            end
            
            [nodes,sdir]=obj.SideDef(sidenum);
            
            J1=obj.Jacobian(pts,elemnum);

            dAx=zeros(length(elemnum),1);
            %% not sure what to do here.
            J1=sqrt(dAx.*dAx+dAy.*dAy);
            %           dAx=sign(sdir).*dAx;
            %           dAy=sign(sdir).*dAy;
        end
        %%%%
    end
    
    methods (Static)
        function [dNx,detJ]=calcGlobalDerivatives(dNx,detJ)
            dNx=dNx./detJ(:,ones(size(dNx,2),1));
        end
    end
end
