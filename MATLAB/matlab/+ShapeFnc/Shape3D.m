classdef Shape3D
    properties
        Connectivity  % nelem x nodes_per_element
        Coordinates   % needs to be consistent with Connectivity numbering so might have all of the nodes
        Dimension = 3
        NumSides 
        NumIntPts 
        NumMassIntPts
        NumSurfIntPts
        MaxIter = 3
        MaxError = 1e-8
    end
    
    properties (Abstract)
        PltConn
    end
    
    methods
        function obj=Shape3D(conn,nodexyz)
            if nargin>0,
                obj.Connectivity=conn;
                switch size(nodexyz,2)
                    case 3
                        obj.Coordinates=nodexyz;
                    case 2
                        obj.Coordinates=[nodexyz zeros(size(nodexyz,1),1)];
                    case 1
                        obj.Coordinates=[nodexyz zeros(size(nodexyz,1),2)];
                end
                    
            end
        end
        function [pts]=Global2Local(obj,PTS,elemnum)
            %pts=obj.LocalCentroid;
            %pts=ones(size(PTS,1),1)*pts;
            pts=obj.initial_guess(PTS,elemnum);
            %pts=repmat(obj.LocalCentroid,length(elemnum),1);
            switch obj.Dimension,
                case 3
                    xiterp=1e16*ones(size(PTS,1),3);
                case 2
                    xiterp=[1e16*ones(size(PTS,1),2) zeros(size(PTS,1),1)];
                    % force zeros in the z's
                    pts(:,3)=zeros(size(pts,1),1);
                    PTS(:,3)=zeros(size(PTS,1),1);
                case 1
                    xiterp=[1e16*ones(size(PTS,1),1) zeros(size(PTS,1),2)];
                    % force zeros in the z's
                    pts(:,2:3)=zeros(size(pts,1),2);
                    PTS(:,2:3)=zeros(size(PTS,1),2);
            end
            err=1e16*ones(size(PTS,1),1);
            niter=1;
            
            Xc=obj.Coordinates(:,1);
            Yc=obj.Coordinates(:,2);
            Zc=obj.Coordinates(:,3);
            
            while (niter<=obj.MaxIter && any(err>obj.MaxError))
                N=obj.Interpolate(pts);
                if length(elemnum)>1,
                    f=PTS-[sum(N.*Xc(obj.Connectivity(elemnum,:)),2) sum(N.*Yc(obj.Connectivity(elemnum,:)),2) sum(N.*Zc(obj.Connectivity(elemnum,:)),2)];
                else
                    f=PTS-[sum(N.*Xc(obj.Connectivity(elemnum,:))',2) sum(N.*Yc(obj.Connectivity(elemnum,:))',2) sum(N.*Zc(obj.Connectivity(elemnum,:))',2)];
                end
                J1=obj.Jacobian(pts,elemnum);
                %% actually need the transpose of J1
                %[ju,det]=VMath.Solve33(J1(:,[1 4 7 2 5 8 3 6 9]),f,obj.Dimension); 
                pts=pts+VMath.Solve33(J1(:,[1 4 7 2 5 8 3 6 9]),f,obj.Dimension); 
                
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
                %error('Maximum iterations exceeded')
            end
        end
        function [PTS]=Local2Global(obj,pts,elemnum)
            if nargin<3,
                elemnum=1:size(obj.Connectivity,1);
            end
            Xc=obj.Coordinates(:,1);
            if obj.Dimension>1,
                Yc=obj.Coordinates(:,2);
            end
            if obj.Dimension>2,
                Zc=obj.Coordinates(:,3);
            end
            N=obj.Interpolate(pts);
            if length(elemnum)>1,
                switch obj.Dimension,
                    case 1
                        PTS=sum(N.*Xc(obj.Connectivity(elemnum,:)),2);
                    case 2
                        PTS=[sum(N.*Xc(obj.Connectivity(elemnum,:)),2) ...
                            sum(N.*Yc(obj.Connectivity(elemnum,:)),2)];
                    case 3
                        PTS=[sum(N.*Xc(obj.Connectivity(elemnum,:)),2) ...
                            sum(N.*Yc(obj.Connectivity(elemnum,:)),2) ...
                            sum(N.*Zc(obj.Connectivity(elemnum,:)),2)];
                    otherwise
                        error('Unknown number of Dimensions: Dim = %d\n',obj.Dimension)
                end
            else
                switch obj.Dimension,
                    case 1
                        PTS=N*Xc(obj.Connectivity(elemnum,:));
                    case 2
                        PTS=[N*Xc(obj.Connectivity(elemnum,:)) ...
                            N*Yc(obj.Connectivity(elemnum,:))];
                    case 3
                        PTS=[N*Xc(obj.Connectivity(elemnum,:)) ...
                            N*Yc(obj.Connectivity(elemnum,:)) ...
                            N*Zc(obj.Connectivity(elemnum,:))];
                    otherwise
                        error('Unknown number of Dimensions: Dim = %d\n',obj.Dimension)
                end
            end
        end
        function [PTS]=GlobalCentroid(obj,elemnum)
            if nargin==1,
                elemnum=1:size(obj.Connectivity,1);
            end
            [pts]=obj.LocalCentroid;
            PTS=obj.Local2Global(repmat(pts,length(elemnum),1),elemnum);
        end
        %%%%
        function [invJ,detJ]=invJacobian(obj,pts,elemnum)
            J1=obj.Jacobian(pts,elemnum);
            [invJ,detJ]=VMath.Inv33(J1,obj.Dimension);
        end
        function [J1,dNx,dNy,dNz]=Jacobian(obj,pts,elemnum)
            switch obj.Dimension
                case 2
                    pts(:,3)=zeros(size(pts,1),1);
                case 1
                    pts(:,2:3)=zeros(size(pts,1),2);
            end
            dNx=obj.lderivX(pts);
            dNy=obj.lderivY(pts);
            dNz=obj.lderivZ(pts);
            x=reshape(obj.Coordinates(obj.Connectivity(elemnum,:),1),length(elemnum),size(obj.Connectivity,2));
            
            switch obj.Dimension
                case 3,
                    y=reshape(obj.Coordinates(obj.Connectivity(elemnum,:),2),length(elemnum),size(obj.Connectivity,2));
                    z=reshape(obj.Coordinates(obj.Connectivity(elemnum,:),3),length(elemnum),size(obj.Connectivity,2));
                case 2
                    y=reshape(obj.Coordinates(obj.Connectivity(elemnum,:),2),length(elemnum),size(obj.Connectivity,2));
                    z=zeros(size(pts,1),length(elemnum));
                case 1
                    y=zeros(size(pts,1),length(elemnum));
                    z=zeros(size(pts,1),length(elemnum));
            end
            
            J1=[sum(x.*dNx,2) sum(y.*dNx,2) sum(z.*dNx,2) ...
                sum(x.*dNy,2) sum(y.*dNy,2) sum(z.*dNy,2)...
                sum(x.*dNz,2) sum(y.*dNz,2) sum(z.*dNz,2)];
            
        end
        function [dNx,dNy,dNz,detJ]=Derivatives(obj,pts,elemnum)
            [J1,dNx,dNy,dNz]=obj.Jacobian(pts,elemnum);
            [dNx,dNy,dNz,detJ]=ShapeFnc.Shape3D.calcGlobalDerivatives(dNx,dNy,dNz,J1,obj.Dimension);
        end
        function [detJ,dAx,dAy,dAz]=SurfJacobian(obj,pts,sidenum,elemnum)
            if nargin==3,
                elemnum=1:size(obj.Connectivity,1);
            end
            
            [nodes,sdir]=obj.SideDef(sidenum);

            J1=obj.Jacobian(pts,elemnum);
            
            dAx=zeros(length(elemnum),1);
            dAy=zeros(length(elemnum),1);
            dAz=zeros(length(elemnum),1);

            switch abs(sdir),
                case 1
                    %dAx=J1(:,5).*J1(:,9)-J1(:,8).*J1(:,6);
                    %dAy=J1(:,8).*J1(:,3)-J1(:,2).*J1(:,9);
                    %dAz=J1(:,2).*J1(:,6)-J1(:,5).*J1(:,3);
                    dAx=J1(:,5).*J1(:,9)-J1(:,6).*J1(:,8);
                    dAy=J1(:,4).*J1(:,9)-J1(:,6).*J1(:,7);
                    dAz=J1(:,4).*J1(:,8)-J1(:,5).*J1(:,7);
                case 2
                    %dAx=J1(:,4).*J1(:,9)-J1(:,7).*J1(:,6);
                    %dAy=J1(:,7).*J1(:,3)-J1(:,1).*J1(:,9);
                    %dAz=J1(:,1).*J1(:,6)-J1(:,4).*J1(:,3);
                    dAx=J1(:,2).*J1(:,9)-J1(:,3).*J1(:,8);
                    dAy=J1(:,1).*J1(:,9)-J1(:,3).*J1(:,7);
                    dAz=J1(:,1).*J1(:,8)-J1(:,2).*J1(:,7);
                case 3
                    %dAx=J1(:,4).*J1(:,8)-J1(:,7).*J1(:,5);
                    %dAy=J1(:,7).*J1(:,2)-J1(:,1).*J1(:,8);
                    %dAz=J1(:,1).*J1(:,5)-J1(:,4).*J1(:,2);
                    dAx=J1(:,2).*J1(:,6)-J1(:,3).*J1(:,5);
                    dAy=J1(:,1).*J1(:,6)-J1(:,3).*J1(:,4);
                    dAz=J1(:,1).*J1(:,5)-J1(:,2).*J1(:,4);
                case 4
                    %% should only be used by wedges
                    Jt=J1(:,1:3)-J1(:,4:6);
                    dAx=Jt(:,2).*J1(:,9) - Jt(:,3).*J1(:,8);
                    dAy=Jt(:,3).*J1(:,7) - Jt(:,1).*J1(:,9);
                    dAz=Jt(:,1).*J1(:,8) - Jt(:,2).*J1(:,7);
                otherwise
                    error('Direction unknown %d in Shape3D.SurfJacobian',abs(sdir));
            end
            detJ=sqrt(VMath.Dot([dAx dAy dAz],[dAx dAy dAz]));
            %           dAx=sign(sdir).*dAx;
            %           dAy=sign(sdir).*dAy;
            %           dAz=sign(sdir).*dAz;
        end
        function [detJ]=calcDetJ(obj,pts,elemnum)
            J1=obj.Jacobian(pts,elemnum);
            detJ=VMath.Det33(J1,obj.Dimension);
        end
        function [vol]=Volume(obj,elemnum)
            if nargin==1,
                elemnum=1:size(obj.Connectivity,1);
            end
            nelem=length(elemnum);
            
            [pts,weights]=obj.getIntPts;
            
            elem=repmat(elemnum,size(pts,1),1);
            
            detJ=reshape(obj.calcDetJ(repmat(pts,nelem,1),elem(:)),size(pts,1),nelem)';
            vol=(detJ*weights);
        end
        
        function [area,elemcent]=SurfArea(obj,sidenums,elemnum)
            area=zeros(length(elemnum),1);
            for i=1:length(elemnum),
                [intpts,weights]=obj.getSurfIntPts(sidenums(i));
                n=length(weights);
                [nodes,sdir]=obj.SideDef(sidenums(i));
                % integrate int(N'*N*dA)
                [detJ]=obj.SurfJacobian(intpts,double(sidenums(i)),elemnum(i)*ones(n,1));
                area(i)=sum(weights.*detJ);
            end
            if nargout==2,
                for i=1:length(elemnum),
                    [cpts]=obj.getSurfIntPts(sidenums(i),1);
                    elemcent(i,:)=obj.Local2Global(repmat(cpts,length(elemnum(i)),1),elemnum(i));
                end
            end
        end
        %%%
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
            a=zeros(nelem,3);
            A=zeros(nelem,1);
            x=obj.Coordinates(:,1);
            y=obj.Coordinates(:,2);
            z=obj.Coordinates(:,3);
            %%
            for i=1:obj.NumSides,
                [pts,we]=obj.getSurfIntPts(i,nintpts);
                elem=repmat(elemnum(:)',size(pts,1),1);
                elem=elem(:);
                sides=obj.SideDef(i);
                sidnod=obj.Connectivity(elemnum,sides);
                N=obj.Interpolate(pts);
                J1=obj.SurfJacobian(repmat(pts,nelem,1),i,elem);
                
                J1r=reshape(repmat(we,nelem,1).*J1,size(pts,1),nelem)';
                %%%%%%
                
                a(:,1)=a(:,1)+sum(J1r.*(x(sidnod)*N(:,sides)'),2);
                a(:,2)=a(:,2)+sum(J1r.*(y(sidnod)*N(:,sides)'),2);
                a(:,3)=a(:,3)+sum(J1r.*(z(sidnod)*N(:,sides)'),2);
                A=A+sum(J1r,2);
            end
            a=a./A(:,ones(1,3));
            
            xx=zeros(nelem,6);
            
            for i=1:obj.NumSides,
                [pts,we]=obj.getSurfIntPts(i,nintpts);
                elem=repmat(elemnum(:)',size(pts,1),1);
                elem=elem(:);
                sides=obj.SideDef(i);
                ns=length(sides);
                sidnod=obj.Connectivity(elemnum,sides);
                N=obj.Interpolate(pts);
                J1=obj.SurfJacobian(repmat(pts,nelem,1),i,elem);
                J1r=reshape(repmat(we,nelem,1).*J1,size(pts,1),nelem)';
                %%%%%%
                xx(:,1)=xx(:,1)+sum(J1r.*(((x(sidnod)-a(:,ones(ns,1)))*N(:,sides)').*((x(sidnod)-a(:,ones(ns,1)))*N(:,sides)')),2);
                xx(:,2)=xx(:,2)+sum(J1r.*(((y(sidnod)-a(:,2*ones(ns,1)))*N(:,sides)').*((y(sidnod)-a(:,2*ones(ns,1)))*N(:,sides)')),2);
                xx(:,3)=xx(:,3)+sum(J1r.*(((z(sidnod)-a(:,3*ones(ns,1)))*N(:,sides)').*((z(sidnod)-a(:,3*ones(ns,1)))*N(:,sides)')),2);
                xx(:,4)=xx(:,4)+sum(J1r.*(((y(sidnod)-a(:,2*ones(ns,1)))*N(:,sides)').*((z(sidnod)-a(:,3*ones(ns,1)))*N(:,sides)')),2);
                xx(:,5)=xx(:,5)+sum(J1r.*(((x(sidnod)-a(:,ones(ns,1)))*N(:,sides)').*((z(sidnod)-a(:,3*ones(ns,1)))*N(:,sides)')),2);
                xx(:,6)=xx(:,6)+sum(J1r.*(((x(sidnod)-a(:,ones(ns,1)))*N(:,sides)').*((y(sidnod)-a(:,2*ones(ns,1)))*N(:,sides)')),2);
            end
            
            B=VMath.SInv33(xx).*A(:,ones(1,6));
            
            if nargout>2,
                nintpts=obj.NumIntPts;
                elem=repmat(elemnum(:)',nintpts*nintpts*nintpts,1);
                elem=elem(:);

                [pts,we]=obj.getIntPts;
                pt=repmat(pts,nelem,1);
                N=obj.Interpolate(pts);

                J1=reshape(obj.calcDetJ(pt,elem),size(pts,1),nelem)';
                Y=zeros(nelem,3);
                Y(:,1)=((x(obj.Connectivity(elemnum,:))*N').*J1)*we;
                Y(:,2)=((y(obj.Connectivity(elemnum,:))*N').*J1)*we;
                Y(:,3)=((z(obj.Connectivity(elemnum,:))*N').*J1)*we;
            end
        end
        %%%
        function plot(obj,elemnum,color)
            if nargin==1 || isempty(elemnum),
                elemnum=1:size(obj.Connectivity,1);
            end
            if nargin<3,
                color='b';
            end
            sz_coord=size(obj.Coordinates,2);
            switch sz_coord
                case 3,
                    plot3(obj.Coordinates(obj.Connectivity(elemnum(1),1),1), ...
                        obj.Coordinates(obj.Connectivity(elemnum(1),1),2), ...
                        obj.Coordinates(obj.Connectivity(elemnum(1),1),3))
                case 2
                    plot(obj.Coordinates(obj.Connectivity(elemnum(1),1),1), ...
                        obj.Coordinates(obj.Connectivity(elemnum(1),1),2))
                case 1
                    plot(obj.Coordinates(obj.Connectivity(elemnum(1),1),1), ...
                        zeros(size(obj.Connectivity(elemnum(1)))))
            end
            vi=obj.PltConn;
            if length(elemnum)>3e3,
                idx=unique(obj.Connectivity(elemnum,:));
                switch sz_coord
                    case 3
                        line(obj.Coordinates(idx,1),obj.Coordinates(idx,2),obj.Coordinates(idx,3), ...
                            'LineStyle','none','Marker','.','Color',color);
                    case 2
                        line(obj.Coordinates(idx,1),obj.Coordinates(idx,2), ...
                            'LineStyle','none','Marker','.','Color',color);
                    case 1
                        line(obj.Coordinates(idx,1),zeros(size(idx)), ...
                            'LineStyle','none','Marker','.','Color',color);
                end
            else
                for i=1:length(elemnum),
                    switch sz_coord
                        case 3
                            line(obj.Coordinates(obj.Connectivity(elemnum(i),vi),1), ...
                                obj.Coordinates(obj.Connectivity(elemnum(i),vi),2), ...
                                obj.Coordinates(obj.Connectivity(elemnum(i),vi),3), ...
                                'Color',color);
                        case 2
                            line(obj.Coordinates(obj.Connectivity(elemnum(i),vi),1), ...
                                obj.Coordinates(obj.Connectivity(elemnum(i),vi),2), ...
                                'Color',color);
                        case 1
                            line(obj.Coordinates(obj.Connectivity(elemnum(i),vi),1), ...
                                zeros(size(obj.Connectivity(elemnum(i),vi))), ...
                                'Color',color);
                    end
                end
            end
        end
        function plot2(obj,elemnum,color)
            if nargin==1,
                elemnum=1:size(obj.Connectivity,1);
            end
            if nargin<3,
                color='b';
            end
            
            nelem=length(elemnum);
            sz_coord=size(obj.Coordinates,2);
 
            x=obj.Coordinates(:,1);
            x=[x(obj.Connectivity(elemnum,:)) NaN(nelem,1)]';
            
            if sz_coord>1,
                y=obj.Coordinates(:,2);
                y=[y(obj.Connectivity(elemnum,:)) NaN(nelem,1)]';
            else
                y=zeros(size(x));
            end
            if sz_coord>2,
                z=obj.Coordinates(:,3);
                z=[z(obj.Connectivity(elemnum,:)) NaN(nelem,1)]';
            else
                z=zeros(size(x));
            end
            %a=[x(:) y(:) z(:)];
            %size(a)
            %length(find(isnan(a(:,1))))
            %a(1:50,:)
            plot3(x(:),y(:),z(:),'Color',color)
        end
        %%
        function [pts]=initial_guess(obj,PTS,elemnum)
            Xc=obj.Coordinates(:,1);
            Yc=obj.Coordinates(:,2);
            Zc=obj.Coordinates(:,3);
            n=size(obj.Connectivity,2);
            
            cpts=obj.GlobalCentroid(elemnum);
           
            % Orgainize points
            if length(elemnum)>1,
                errx=[repmat(PTS(:,1),1,n)-Xc(obj.Connectivity(elemnum,:)) cpts(:,1)];
                erry=[repmat(PTS(:,2),1,n)-Yc(obj.Connectivity(elemnum,:)) cpts(:,2)];
                errz=[repmat(PTS(:,3),1,n)-Zc(obj.Connectivity(elemnum,:)) cpts(:,3)];
            else
                errx=[repmat(PTS(:,1),1,n)-Xc(obj.Connectivity(elemnum,:))' cpts(:,1)];
                erry=[repmat(PTS(:,2),1,n)-Yc(obj.Connectivity(elemnum,:))' cpts(:,2)];
                errz=[repmat(PTS(:,3),1,n)-Zc(obj.Connectivity(elemnum,:))' cpts(:,3)];
            end
            [ju,I]=min(errx.*errx+erry.*erry+errz.*errz,[],2);
            hpts=[obj.LocalNode;obj.LocalCentroid];
            pts=hpts(I,:);
        end
        
        %%%
        function [uelems,usides,unperside]=UniqueSideNum(obj)
            nmin=min(obj.Connectivity(:))-1;
            nnod=length(unique(obj.Connectivity(:)));
            n=size(obj.Connectivity,1);
            sides=zeros(obj.NumSides*n,1);
            elems=zeros(obj.NumSides*n,1);
            unm=zeros(obj.NumSides*n,1);
            nnodesperside=zeros(obj.NumSides*n,1);
            nn=zeros(n,1);
            
            for i=1:obj.NumSides,
                [sid]=obj.SideDef(i);

                node=double(obj.Connectivity(:,sid)-nmin);
                
                [nm,idx]=min(node,[],2);
                
                idx=mod(idx+2,length(sid));
               
                if any(idx==0),
                    idx(idx==0)=length(sid);
                end
                for j=1:n,
                    nn(j)=node(j,idx(j));
                end
                
                unm(((i-1)*n+1):(i*n))=nn+nm*nnod;
               
                sides(((i-1)*n+1):(i*n))=i*ones(n,1);
                elems(((i-1)*n+1):(i*n))=(1:n)';
                nnodesperside(((i-1)*n+1):(i*n))=length(sid)*ones(n,1);
            end
            [unm,idx1]=sort(unm);
            sides=sides(idx1);
            elems=elems(idx1);
            nnodesperside=nnodesperside(idx1);
            [unum,idx]=setdiff(unm,unique(unm(~diff(unm))));
            usides=sides(idx);
            uelems=elems(idx);
            [uelems usides unum];
            unperside=nnodesperside(idx);
        end
    end
    
    methods (Abstract, Static)
        N=Interpolate(pts)
        [pts]=LocalCentroid
        [xb,yb,zb]=LocalBounds
        dNx=lderivX(pts)
        dNy=lderivY(pts)
        dNz=lderivZ(pts)
        [node,sdir]=SideDef(sidenum)
        [v,conn]=GenSubCells(npts,type)
    end
    
    methods (Abstract)
        [pts,weights]=getIntPts(obj,n)
        [pts,weights]=getSurfIntPts(obj,side,n)
        [pts,weights]=getGenIntPts(obj,n,limits)
    end
    
    methods (Static)
        function [intpts,weight,idnums]=GaussWeights(npts)
            %
            %  gaussw gives the local coordinates and the cooresponding weights
            %	for use in the shape function routines
            %
            %
            %	function [intpts,weight,idnums]=gaussw(npts)
            %
            %   [intpts,weight,idnums]=gaussw([3 2 1]);
            %
            %
            ndim=length(npts);
            pts=zeros(max(npts),ndim);
            idnum=zeros(max(npts),ndim);
            we=pts;
            for i=1:ndim,
                if npts(i)==1,
                    pts(1,i)=0;
                    we(1,i)=2;
                elseif npts(i)==2,
                    a=0.57735026918963;
                    pts(1:2,i)=[-a;a];
                    idnum(1:2,i)=[0;1];
                    we(1:2,i)=[1;1];
                elseif npts(i)==3,
                    a=0.77459666924148;
                    pts(1:3,i)=[-a;0;a];
                    idnum(1:3,i)=[0;1;2];
                    a=0.555555555555555;
                    b=0.88888888888889;
                    we(1:3,i)=[a;b;a];
                elseif npts(i)==4,
                    a=0.86113631159495;
                    b=0.33998104358486;
                    pts(1:4,i)=[-a;-b;b;a];
                    idnum(1:4,i)=[0;1;2;3];
                    a=0.347854845137454;
                    b=0.652145154862545;
                    we(1:4,i)=[a;b;b;a];
                elseif npts(i)==5,
                    a=0.90617984593866;
                    b=0.53846931010568;
                    pts(1:5,i)=[-a;-b;0;b;a];
                    idnum(1:5,i)=[0;1;2;3;4];
                    a=0.23692688505619;
                    b=0.478628670499366;
                    c=0.568888888888889;
                    we(1:5,i)=[a;b;c;b;a];
                elseif npts(i)==6,
                    a=0.93246951420315;
                    b=0.66120938646627;
                    c=0.23861918608320;
                    pts(1:6,i)=[-a;-b;-c;c;b;a];
                    idnum(1:6,i)=[0;1;2;3;4;5];
                    a=0.17132449237917;
                    b=0.36076157304814;
                    c=0.467913934572691;
                    we(1:6,i)=[a;b;c;c;b;a];
                elseif npts(i)==7,
                    a=0.94910791234276;
                    b=0.74153118559939;
                    c=0.40584515137740;
                    pts(1:7,i)=[-a;-b;-c;0;c;b;a];
                    idnum(1:7,i)=[0;1;2;3;4;5;6];
                    a=0.12948496616887;
                    b=0.279705391489277;
                    c=0.381830050505119;
                    d=0.417959183673469;
                    we(1:7,i)=[a;b;c;d;c;b;a];
                elseif npts(i)==8,
                    a=0.96028985649754;
                    b=0.79666647741363;
                    c=0.52553240991633;
                    d=0.18343464249565;
                    pts(1:8,i)=[-a;-b;-c;-d;d;c;b;a];
                    idnum(1:8,i)=[0;1;2;3;4;5;6;7];
                    a=0.10122853629038;
                    b=0.22238103445337;
                    c=0.31370664587789;
                    d=0.36268378337836;
                    we(1:8,i)=[a;b;c;d;d;c;b;a];
                elseif npts(i)==9,
                    a=0.96816023950763;
                    b=0.83603110732664;
                    c=0.61337143270059;
                    d=0.32425342340381;
                    pts(1:9,i)=[-a;-b;-c;-d;0;d;c;b;a];
                    idnum(1:9,i)=[0;1;2;3;4;5;6;7;8];
                    a=0.08127438836157;
                    b=0.18064816069486;
                    c=0.26061069640294;
                    d=0.31234707704000;
                    e=0.33023935500126;
                    we(1:9,i)=[a;b;c;d;e;d;c;b;a];
                elseif npts(i)==10,
                    a=0.97390652851717;
                    b=0.86506336668899;
                    c=0.67940956829902;
                    d=0.43339539412925;
                    e=0.14887433898163;
                    pts(1:10,i)=[-a;-b;-c;-d;-e;e;d;c;b;a];
                    idnum(1:10,i)=[0;1;2;3;4;5;6;7;8;9];
                    a=0.066671344308688;
                    b=0.149451349150581;
                    c=0.219086362515982;
                    d=0.269266719309996;
                    e=0.295524224714753;
                    we(1:10,i)=[a;b;c;d;e;e;d;c;b;a];
                else
                    u=(1:npts(i)-1)./sqrt((2*(1:npts(i)-1)).^2-1);
                    [vc,gp]=eig(diag(u,-1)+diag(u,1));
                    [pts(1:npts(i),i),k]=sort(diag(gp));
                    idnum=[];
                    we(1:npts(i),i)=2*vc(1,k)'.^2;
                end
            end
            switch ndim
                case 1
                    intpts=pts;
                    weight=conj(we);
                    idnums=idnum;
                case 2
                    [px,py]=meshgrid(pts(1:npts(1),1),pts(1:npts(2),2));
                    [wts1,wts2]=meshgrid(we(1:npts(1),1),we(1:npts(2),2));
                    if ~isempty(idnum),
                        [idn1,idn2]=meshgrid(idnum(1:npts(1),1),idnum(1:npts(2),2));
                        idnums=[idn1(:) idn2(:)];
                    end
                    weight=wts1(:).*wts2(:);
                    intpts=[px(:) py(:)];
                case 3
                    [px,py,pz]=meshgrid(pts(1:npts(1),1),pts(1:npts(2),2),pts(1:npts(3),3));
                    [wts1,wts2,wts3]=meshgrid(we(1:npts(1),1),we(1:npts(2),2),we(1:npts(3),3));
                    if ~isempty(idnum),
                        [idn1,idn2,idn3]=meshgrid(idnum(1:npts(1),1),idnum(1:npts(2),2),idnum(1:npts(3),3));
                        idnums=[idn1(:) idn2(:) idn3(:)];
                    end
                    weight=wts1(:).*wts2(:).*wts3(:);
                    intpts=[px(:) py(:) pz(:)];
                otherwise
                    error('Unknown dimension in GaussWeights')
            end   
        end
        function [Points,Weights,idnums]=TriIntegration(npts)
            switch npts %% npts is actually the order of integration
                case 1   % 1st order
                    Points=[1 1]/3;
                    Weights=1;
                    
                case 2    % second order
                    Points=[1/6 1/6
                        2/3 1/6
                        1/6 2/3];
                    Weights=[1 1 1]'/3;
                case -2   % third order
                    Points=[1/3 1/3;
                        0.6 0.2
                        0.2 0.6
                        0.2 0.2];
                    Weights=[-27 25 25 25]'/48;
                    
                case 4  % fourth order
                    a=0.109951743655322;
                    b=0.223381589678011;
                    a1=0.816847572980459;
                    b1=0.091576213509771;
                    a2=0.108103018168070;
                    b2=0.445948490915965;
                    Points=[a1 b1
                        b1 a1
                        b1 b1
                        a2 b2
                        b2 a2
                        b2 b2];
                    Weights=[a;a;a;b;b;b];
                case 3   % fifth order
                    a=0.225;
                    b=0.132394152788506;
                    c=0.125939180544827;
                    a1=0.059715871789770;
                    b1=0.470142064105115;
                    a2=0.797426985353087;
                    b2=0.101286507323456;
                    Points=[1/3 1/3;
                        a1 b1
                        b1 a1
                        b1 b1
                        a2 b2
                        b2 a2
                        b2 b2];
                    Weights=[a;b;b;b;c;c;c];
                case -400
                    %% cludge for Wedge15 elements
                    a1=0.5;
                    Points=[a1 a1
                        0 a1
                        a1 0];
                    Weights=[1 1 1]'/3;
                case -500
                    %% cludge for Wedge15 elements
                    a3=0.81684757298046;
                    a4=0.09157621350977;
                    a5=0.10810301816807;
                    a6=0.44594849091597;
                    Points=[a3 a4
                        a4 a3
                        a4 a4
                        a5 a6
                        a6 a5
                        a6 a6];
                    a1=0.1099517436553321;
                    a2=0.2233815896780282;
                    Weights=[a1 a1 a1 a2 a2 a2]';
                case 5  %% order 6
                    a=0.050844906370207;
                    b=0.116786275726379;
                    c=0.082851075618374;
                    a1=0.873821971016996;
                    b1=0.063089014491502;
                    a2=0.501426509658179;
                    b2=0.249286745170910;
                    a3=0.636502499121399;
                    b3=0.310352451033785;
                    c3=0.053145049844816;
                    Points=[a1 b1
                        b1 a1
                        b1 b1
                        a2 b2
                        b2 a2
                        b2 b2
                        a3 b3
                        a3 c3
                        b3 a3
                        b3 c3
                        c3 a3
                        c3 b3];
                    Weights=[a;a;a;b;b;b;c;c;c;c;c;c];
                case 6  % 7th order
                    a=-0.149570044467670;
                    b=0.175615257433204;
                    c=0.053347235608839;
                    d=0.077113760890257;
                    a2=0.479308067841923;
                    b2=0.260345966079038;
                    a3=0.869739794195568;
                    b3=0.065130102902216;
                    a4=0.638444188569809;
                    b4=0.312865496004875;
                    c4=0.048690315425316;
                    Points=[1/3 1/3
                        a2 b2
                        b2 a2
                        b2 b2
                        a3 b3
                        b3 a3
                        b3 b3
                        a4 b4
                        a4 c4
                        b4 a4
                        b4 c4
                        c4 a4
                        c4 b4];
                    Weights=[a;b;b;b;c;c;c;d;d;d;d;d;d];
                otherwise
                    error('Number of Points not implemented')
            end
            Weights=0.5*Weights;  %% this is necessary because of how Cowper defined the weights.
            %  He defined them relative to the Area while the rest of the
            %  code here is dependent on the jacobian which is twice the
            %  area in a triangle.
            idnums=ShapeFnc.Shape3D.genGaussPointNames(Points);
        end
        function [Points,Weights,idnums]=TetIntegration(npts)
            switch npts
                case 1 %% order 1
                    Points=[1 1 1]/4;
                    Weights=1;
                case 2  %% order 2
                    a=0.585410196624969;
                    b=0.13819660112501;
                    Points=[a b b
                        b a b
                        b b a
                        b b b];
                    Weights=ones(4,1)/4;
                    
                case 3  %% Order 3
                    a=0.0948472649145;
                    b=0.24127699682327;
                    c=0.569028473347700;
                    Points=[a a b
                        a a c
                        a b a
                        a b c
                        a c a
                        a c b
                        b a a
                        b a c
                        b c a
                        c a a
                        c a b
                        c b a];
                    Weights=ones(12,1)/12;
                    
                case 4
                    a=0.04254602077708200;
                    b=0.1126879257180180;
                    c=0.07349304311636401;
                    a1=0.454496295874350;
                    b1=0.045503704125650;
                    a2=0.310885919263301;
                    b2=0.067342242210098;
                    a3=0.721794249067326;
                    b3=0.092735250310891;
                    Points=[a1 a1 b1
                        a1 b1 a1
                        a1 b1 b1
                        b1 a1 a1
                        b1 a1 b1
                        b1 b1 a1
                        a2 a2 a2
                        b2 a2 a2
                        a2 b2 a2
                        a2 a2 b2
                        b3 b3 b3
                        a3 b3 b3
                        b3 a3 b3
                        b3 b3 a3];
                    Weights=[a;a;a;a;a;a;b;b;b;b;c;c;c;c];
                    
                case 5
                    a=0.00839563235;
                    b=0.01109034477167;
                    a1=.77164290200000;
                    b1=.07611903264000;
                    a2=.40423391340000;
                    b2=0.1197005277000;
                    c2=0.07183164526000;
                    Points=[a1 b1 b1
                        b1 a1 b1
                        b1 b1 a1
                        %
                        b1 b1 b1
                        a2 c2 b2
                        a2 b2 c2
                        %
                        b2 c2 a2
                        b2 a2 c2
                        c2 b2 a2
                        %
                        c2 a2 b2
                        a2 c2 a2
                        c2 a2 a2
                        %
                        a2 a2 c2
                        b2 a2 a2
                        a2 b2 a2
                        a2 a2 b2];
                        Weights=[a;a;a;a;b;b;b;b;b;b;b;b;b;b;b;b];
                        
                otherwise
                    error('Number of Points not implemented')
            end
            Weights=Weights/6;
            idnums=ShapeFnc.Shape3D.genGaussPointNames(Points);
        end
        
        function [Points,Weights,idnums]=WedIntegration(npts)
            switch npts,
                case 1
                    Points=[1/3 1/3 0];
                    Weights=1;
                case 2
                    a=1/6;
                    b=1/sqrt(3);
                    c=2/3;
                    Points=[a a -b
                        c a -b
                        a c -b
                        a a b
                        c a b
                        a c b];
                    Weights=ones(6,1)/6;
                case 4
                    a1=0.5;
                    a2=0.77459666924148;
                    Points=[a1 a1 -a2
                        0 a1 -a2
                        a1 0 -a2
                        a1 a1 0
                        0 a1 0
                        a1 0 0
                        a1 a1 a2
                        0 a1 a2
                        a1 0 a2];
                    a1=0.09259259259259;
                    a2=0.14814814814815;
                    Weights=[a1 a1 a1 a2 a2 a2 a1 a1 a1]';
                case 5
                    a2=0.77459666924148;
                    a3=0.81684757298046;
                    a4=0.09157621350977;
                    a5=0.10810301816807;
                    a6=0.44594849091597;
                    Points=[a3 a4 -a2
                        a4 a3 -a2
                        a4 a4 -a2
                        a5 a6 -a2
                        a6 a5 -a2
                        a6 a6 -a2
                        a3 a4 0
                        a4 a3 0
                        a4 a4 0
                        a5 a6 0
                        a6 a5 0
                        a6 a6 0
                        a3 a4 a2
                        a4 a3 a2
                        a4 a4 a2
                        a5 a6 a2
                        a6 a5 a2
                        a6 a6 a2];
                    a1=0.03054215101537;
                    a2=0.04886744162459;
                    a3=0.06205044157723;
                    a4=0.09928070652356;
                    Weights=[a1 a1 a1 a3 a3 a3 a2 a2 a2 a4 a4 a4 a1 a1 a1 a3 a3 a3]';
                otherwise
                    error('Uknown number of gauss points')
            end
            idnums=ShapeFnc.Shape3D.genGaussPointNames(Points);
        end
        function [dNx,dNy,dNz,detJ]=calcGlobalDerivatives(dNx,dNy,dNz,J1,dim)
            f=zeros(size(dNx,1),3,size(dNx,2));
            f(:,1,:)=dNx;
            f(:,2,:)=dNy;
            f(:,3,:)=dNz;
            
            [sol,detJ]=VMath.Solve33(J1,f,dim);
            dNx=squeeze(sol(:,1,:));
            dNy=squeeze(sol(:,2,:));
            dNz=squeeze(sol(:,3,:));
                
        end      
        function [intpts,weight]=GeneralGaussWeights(npts,limits)
            %
            %   limits=[minx maxx;miny maxy;minz maxz];
            %
            [intpts,weight]=ShapeFnc.Shape3D.GaussWeights(npts);
            intpts(:,1)=intpts(:,1)*(limits(1,2)-limits(1,1))/2+(limits(1,1)+limits(1,2))/2;
            weight=weight*(limits(1,2)-limits(1,1))/2;
            if size(intpts,2)>1,
                intpts(:,2)=intpts(:,2)*(limits(2,2)-limits(2,1))/2+(limits(2,1)+limits(2,2))/2;
                weight=weight*(limits(2,2)-limits(2,1))/2;
            end
            if size(intpts,2)>2,
                intpts(:,3)=intpts(:,3)*(limits(3,2)-limits(3,1))/2+(limits(3,1)+limits(3,2))/2;
                weight=weight*(limits(3,2)-limits(3,1))/2;
            end
        end
        function [gp_id]=genGaussPointNames(gp)
            [n,dim]=size(gp);
            gp_id=zeros(n,dim);
            
            [a1,I1,J1]=unique(gp(:,1));
            gh1=(0:length(a1))';
            gp_id(:,1)=gh1(J1);
            
            if dim>1,
                [a2,I2,J2]=unique(gp(:,2));
                gh2=(0:length(a2))';
                gp_id(:,2)=gh2(J2);
            end
            if dim>2,
                [a3,I3,J3]=unique(gp(:,3));
                gh3=(0:length(a3))';
                gp_id(:,3)=gh3(J3);
            end
            
            
        end
        
    end
end
