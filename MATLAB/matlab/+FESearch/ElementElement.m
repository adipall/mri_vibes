classdef ElementElement
    properties
        ModelObj
        BlocksOfInterest
        Dimension = 3
        Tol = 0.05  % Amount to allow for extrapolation
        SearchTol = []  % use default in the search
        TightTol = 0.05 % amount where it doesn't matter whether you extrapolate or interpolate
    end
    methods
        function obj=ElementElement(fexo1,idxblk1,tol,stol)
            if nargin>0,
                obj.ModelObj=fexo1;
                obj.BlocksOfInterest=idxblk1;
                if nargin>3 &&~isempty(tol)
                    obj.Tol=tol;
                end
                if nargin>4 && ~isempty(stol),
                    obj.SearchTol=stol;
                end
                obj.Dimension=size(fexo1.Nodes.Coordinates,2);
            end
        end
        function [Vol,Ielem2,vol1,vol2]=ElementIntersectionClark(obj,fexo,idxblk)
            
            vol1=[];
            vol1s=[];
            nelem1=1;
            PTS1=[];
            for i=1:length(obj.BlocksOfInterest),
                blkid=obj.BlocksOfInterest(i);
                ne=size(obj.ModelObj.Blocks(blkid).Connectivity,1);
                Sobj1{i}=ShapeFactory.CreateShape(obj.ModelObj.Blocks(blkid).ElementType, ...
                    obj.ModelObj.Blocks(blkid).Connectivity,obj.ModelObj.Nodes.Coordinates);
                pts=Sobj1{i}.getIntPts;
                npts=size(pts,1);
                elem1=repmat((1:ne),npts,1);
                elem1=elem1(:);
                PTSa=Sobj1{i}.Local2Global(repmat(pts,ne,1),elem1);
                volsub=repmat((Sobj1{i}.Volume)'/npts,npts,1);
                vol1s=cat(1,vol1s,volsub(:));
                vola=Sobj1{i}.Volume;
                vol1=cat(1,vol1,vola(:));
                PTS1=cat(1,PTS1,PTSa);
                nelem1=cat(2,nelem1,ne);
            end
            vol2=[];
            
            nelem2=1;
            for i=1:length(idxblk),
                blkid=idxblk(i);
                Sobj2{i}=ShapeFactory.CreateShape(fexo.Blocks(blkid).ElementType, ...
                    fexo.Blocks(blkid).Connectivity,fexo.Nodes.Coordinates);
                
                vola=Sobj2{i}.Volume;
                vol2=cat(1,vol2,vola);
                nelem2=cat(2,nelem2,size(fexo.Blocks(blkid).Connectivity,1));
            end
            SeaObj=FESearch.ElementNode(fexo,idxblk);
           
            [lcoords,elems,blks,nodes,nodesnotfound]=SeaObj.LocalCoord(PTS1,scale);
            
            %%%%%%%%%%%
            Vol=cell(sum(nelem2)-1,1);
            Ielem2=cell(sum(nelem2)-1,1);
            for i=1:length(Ielem2),
                idx=find(elems==i)
                if ~isempty(idx),
                    el=elem1(idx);
                    [ho,ind]=unique(el);
                    vv=zeros(length(ind),1);
                    for j=1:length(ind),
                        vv(j)=sum(vol1s(nodes(ind(j)==ind)));
                    end
                    Ielem2{i}=ho;
                    Vol{i}=vv*vol2(i)/sum(vv);
                end
                %keyboard
            end
            %keyboard
            for i=1:length(Vol),
                Vol{i}=Vol{i};
            end
            vol1=vol1;
            vol2=vol2;
        end
        
        function [Vol,Ielem2,vol1,vol2]=ElementIntersectionRashid(obj,fexo,idxblk)
            
            %dim=size(fexo.Nodes.Coordinates,2);
            B1=[];
            a1=[];
            vol1=[];
            nelem1=1;
            for i=1:length(obj.BlocksOfInterest),
                blkid=obj.BlocksOfInterest(i);
                
                Sobj1{i}=ShapeFactory.CreateShape(obj.ModelObj.Blocks(blkid).ElementType, ...
                    obj.ModelObj.Blocks(blkid).Connectivity,obj.ModelObj.Nodes.Coordinates);
                [Ba,aa]=Sobj1{i}.QuadApprox;
                vola=Sobj1{i}.Volume;
                B1=cat(1,B1,Ba);
                a1=cat(1,a1,aa);
                vol1=cat(1,vol1,vola);
                nelem1=cat(2,nelem1,size(obj.ModelObj.Blocks(blkid).Connectivity,1));
            end
            B2=[];
            a2=[];
            Y=[];
            vol2=[];
            nelem2=1;
            for i=1:length(idxblk),
                blkid=idxblk(i);
                Sobj{i}=ShapeFactory.CreateShape(fexo.Blocks(blkid).ElementType, ...
                    fexo.Blocks(blkid).Connectivity,fexo.Nodes.Coordinates);
                
                %[Ba,aa,Ya]=Sobj(i).QuadApprox;
                [Ba,aa]=Sobj{i}.QuadApprox;
                vola=Sobj{i}.Volume;
                B2=cat(1,B2,Ba);
                a2=cat(1,a2,aa);
                %Y=cat(1,Y,Ya);
                vol2=cat(1,vol2,vola);
                nelem2=cat(2,nelem2,size(fexo.Blocks(blkid).Connectivity,1));
            end
            %%%%%
            
            SeaObj=FESearch.Search;
            if obj.Dimension==2,
                coord1=[fexo.Nodes.Coordinates zeros(size(fexo.Nodes.Coordinates,1),1)];
                coord2=[obj.ModelObj.Nodes.Coordinates zeros(size(obj.ModelObj.Nodes.Coordinates,1),1)];
            else
                coord1=fexo.Nodes.Coordinates;
                coord2=obj.ModelObj.Nodes.Coordinates;
            end
  %%%%%%%%%%%%%%%%%%%%%%% REMOVE%%%%%%%%%%%%%%%%%%%%%%
            if 1,
                %Ielem1p=SeaObj.BlockElementSearch(coord2, ...
                %    obj.ModelObj.Blocks(obj.BlocksOfInterest),coord1, ...
                %    fexo.Blocks(idxblk), ...
                %    [0 0 0]);
                Ielem2=SeaObj.BlockElementSearch(coord1, ...
                    fexo.Blocks(idxblk),coord2, ...
                    obj.ModelObj.Blocks(obj.BlocksOfInterest), ...
                    [0 0 0]);
                %[Ielem1,Ielem2]=FESearch.ElementElement.setdiff_cell(Ielem1p,Ielem2p);
                %plot(obj.ModelObj,fexo,1,Ielem2)
                %save data2 Ielem2
            else
                %load data
                load data
            end
            %%%%%%%%%%%
            
            Vol=cell(sum(nelem2)-1,1);
            parfor i=1:length(Ielem2),
                %fprintf('%d\n',i)
                if isempty(Y),
                    YY=[];
                else
                    YY=Y(i,:);
                end
                if ~isempty(Ielem2{i}),
                    switch obj.Dimension
                        case 3
                            %[voll,idx]=obj.calcVolRashid3d(B2(i,:),a2(i,:),YY,vol2(i),B1(Ielem2{i},:),a1(Ielem2{i},:),i,Ielem2{i});
                            [voll,idx]=obj.calcVolRashid3d_old(B2(i,:),a2(i,:),YY,vol2(i),B1(Ielem2{i},:),a1(Ielem2{i},:),i,Ielem2{i});
                            %vols=[];idx=Ielem2{i};
                        case 2
                            [voll,idx]=obj.calcVolRashid2d(B2(i,:),a2(i,:),YY,vol2(i),B1(Ielem2{i},:),a1(Ielem2{i},:),i,Ielem2{i});
                    end
                    if any(voll<0) || any(isnan(voll)) || isempty(voll),
                        error('Negative or nan vols')
                    end
                    Ielem2{i}=Ielem2{i}(idx);
                    Vol{i}=voll;
                else
                    Vol{i}=[];
                end
            end
            %keyboard
            for i=1:length(Vol),
                Vol{i}=Vol{i};
            end
            vol1=vol1;
            vol2=vol2;
        end
    end
    
    methods (Static)
        function [vols,vi]=calcVolRashid2d(B,a,Y,vol,B1,a1,i,Ielem2)
            if size(B,1) >1,
                error('calcVolRashid can only determine the volume fractions in one element')
            end
            C=repmat(B,size(B1,1),1)+B1;
            %%
            [Cinv,detC]=VMath.SInv22(C);
            ba=VMath.SAx22(B,a);
            %%
            b1a=VMath.SAx22(B1,a1);
            %%
            c=repmat(ba,size(b1a,1),1)+b1a;
            c=VMath.SAx22(Cinv,c);
            %%%
            ca=VMath.SAx22(C,c);
            bc=repmat(a(1)*ba(1)+a(2)*ba(2),size(B1,1),1) ...
                + (a1(:,1).*b1a(:,1)+a1(:,2).*b1a(:,2)) ...
                - (c(:,1).*ca(:,1)+c(:,2).*ca(:,2));
            %%
            v=c;
            vi=1:size(B1,1);
            
            flag=1;
            while flag,
                vols=FESearch.ElementElement.calcVolRashid2dIter(C,c,bc,detC,v,vol,vi,i,Ielem2);
                if any(vols<0),
                    [ju,idx]=min(vols(vi));
                    vi(idx)=[];
                    vols(idx)=0;
                else
                    flag=0;
                end
            end
            %%% apply the correction
            if ~isempty(Y),
                for i=1:3,
                    vols=FESearch.ElementElement.calcVolRashid2dIter(C,c,bc,detC,v,vol,vi,i,Ielem2);
                    A=sum(repmat(vols(vi,:),1,3).*Cinv(vi,:),1);
                    
                    a=Y-sum(repmat(vols(vi),1,2).*c(vi,:),1);
                    
                    a=VMath.SSolve22(A,a);
                    v=VMath.SAx22(Cinv(vi,:),a) + c;
                end
            end
        end
        
        function [vols,vi]=calcVolRashid3d(B,a,Y,vol,B1,a1,i,Ielem2)
            if size(B,1) >1,
                error('calcVolRashid can only determine the volume fractions in one element')
            end
            C=repmat(B,size(B1,1),1)+B1;
            %%
            %[Cinv,detC]=VMath.SInv33(C);

            ba=VMath.SAx33(B,a);
            %%
            b1a=VMath.SAx33(B1,a1);
            %%
            ca=repmat(ba,size(b1a,1),1)+b1a;
            %c=VMath.SAx33(Cinv,c);
            c=VMath.SSolve33(C,ca);
            %%%
            %ca=VMath.SAx33(C,c);
            bc=repmat(VMath.Dot(a,ba),size(B1,1),1) ...
                + VMath.Dot(a1,b1a) - VMath.Dot(c,ca);
            %
            v=c;
            
            flag=1;
            vi=1:size(B1,1);
            % get an initial state
            %Ul=ones(size(B1,1),1);
            %Ul=normrnd(0,1,size(B1,1),1);
            %Ul=(vol/sum(Ul))*Ul;
            [Ul,ncflag]=FESearch.ElementElement.calcVolRashid3dIter(C,c,bc,detC,v,vol,vi,i,Ielem2);
            if min(real(Ul))<0,
                [ju,idx]=min(real(Ul));
                Ul(idx)=[];
                vi(idx)=[];
            end
            restart=0;
            %vols=[];
            nrestarts=10;
            while flag && restart<nrestarts,
                [W,ncflag]=FESearch.ElementElement.calcVolRashid3dIter(C,c,bc,detC,v,vol,vi,i,Ielem2);
                if length(vi)>1,
                    W=W-Ul;
                    t=ones(length(W),1);
                    idx=find(real(W)<0);
                    
                    t(idx)=-Ul(idx)./W(idx);
                    [mt,id]=min(t);
                    Ul=Ul+mt*W;
                else
                    idx=[];
                    if Ul<0;
                        Ul=0;
                    end
                end
                if ~isempty(idx), %abs(mt-1)>0,
                    vi(id)=[];
                    Ul(id)=[];
                else
                    volshold=real(Ul);
                    if isempty(volshold),
                        Ul
                        Ielem2
                        i
                    end
                    vols=volshold;
                    
                    
                    %if ~isempty(idx),
                    %   disp('negative')
                    %end
                    %if max(abs(vols(idx)))/vol<1e-4,
                    vols(volshold<0)=[];
                    vi(volshold<0)=[];
                    %else
                    %    pause;
                    %end
                    %length(vols)
                    if abs(sum(vols)-vol)/vol >1e-3,
                        %[mt sum(vols) vol]
                    end
                    flag=0;
                    if sum(vols)~=0 && ~isempty(vols),
                        vols=vol/sum(vols)*vols;
                    else
                        Ul
                        vols=[];
                        Ul=[];
                    end
                end
                if (ncflag && restart<nrestarts )|| isempty(real(Ul)),
                    % Try a restart
                    vi=1:size(B1,1);
                    % get an initial state
                    %Ul=ones(size(B1,1),1);
                    Ul=normrnd(0,1,size(B1,1),1);
                    Ul=(vol/sum(Ul))*Ul;
                    ncflag=0;
                    flag=1;
                    restart=restart+1;
                    fprintf('%d ******************\n',i)
                end
            end
            if ncflag || restart==nrestarts,
                W
                error('Newton loops for element %d did not converge',i)
            end
            
            %%% apply the correction
            
            if ~isempty(Y),
                holdvols=vols;
                for i=1:3,
                    vols=FESearch.ElementElement.calcVolRashid3dIter(C,c,bc,detC,v,vol,vi,i,Ielem2);
                    vols=real(vols);
                    if any(vols<0),
                        warning('Negative Volumes')
                    end
                    A=sum(repmat(vols,1,6).*Cinv(vi,:),1);
                    
                    %[Ainv]=VMath.SInv33(A);

                    a=Y-sum(repmat(vols,1,3).*c(vi,:),1);
                    
                    a=VMath.SSolve33(A,a);
                    %v=VMath.SAx33(Cinv,a) + c;
                    v=VMath.SSolve33(C,a)+c;
                end
            end
           
        end
        function [vols,vi]=calcVolRashid3d_old(B,a,Y,vol,B1,a1,i,Ielem2)
            if size(B,1) >1,
                error('calcVolRashid can only determine the volume fractions in one element')
            end
            C=repmat(B,size(B1,1),1)+B1;
            %%
            %[Cinv,detC]=VMath.SInv33(C);

            ba=VMath.SAx33(B,a);
            %%
            b1a=VMath.SAx33(B1,a1);
            %%
            ca=repmat(ba,size(b1a,1),1)+b1a;
            [c,detC]=VMath.SSolve33(C,ca);
            %%%
            %ca=VMath.SAx33(C,c);
            bc=repmat(VMath.Dot(a,ba),size(B1,1),1) ...
                + VMath.Dot(a1,b1a) - VMath.Dot(c,ca);
            
            %
            v=c;
            vi=1:size(B1,1);
            %vi=find(bc>0);
            flag=1;
            %vols=ones(size(B1,1));
            while flag,
                vols=FESearch.ElementElement.calcVolRashid3dIter(C,c,bc,detC,v,vol,vi,i,Ielem2);
                if any(real(vols)<0) %|| any(imag(vols)>1e-12),
                    idx=find(real(vols)<0);
                    vi(idx)=[];
                    vols(idx)=[];
                else
                    if any(real(vols)/max(real(vols))<1e-4),
                        idx=find(real(vols)/max(real(vols))<1e-3);
                        vi(idx)=[];
                        vols(idx)=[];
                    else
                        if abs(max(imag(vols))/max(real(vols)))>1e-5,
                            warning('Imaginary Volumes')
                        end
                        vols=real(vols);
                        flag=0;
                    end
                end
            end
            %%% apply the correction
            
            if ~isempty(Y),
                holdvols=vols;
                for i=1:3,
                    vols=FESearch.ElementElement.calcVolRashid3dIter(C,c,bc,detC,v,vol,vi,i,Ielem2);
                    vols=real(vols);
                    if any(vols<0),
                        warning('Negative Volumes')
                        holdvols
                        vols
                    end
                    A=sum(repmat(vols,1,6).*Cinv(vi,:),1);
                    
                    %[Ainv]=VMath.SInv33(A);

                    a=Y-sum(repmat(vols,1,3).*c(vi,:),1);
                    
                    a=VMath.SSolve33(A,a);
                    v=VMath.SSolve33(C,a) + c;
                end
            end
           
        end
        function [vols]=calcVolRashid2dIter(C,c,bc,detC,v,vol,vi,i,Ielem2)
            
            vols=zeros(size(C,1),1); %% until I figure out better
            s=1;  % until I figure out better

            va=v(vi,:)-c(vi,:);
            
            P=0.5/pi*(detC(vi).^(1/2));
            
            Q=sum(va.*VMath.SAx22(C(vi,:),va),2) + bc(vi);

            lam=(vol +sum(Q./(2*P)))/sum(1./(2*P));
            
            vols(vi)=vols(vi)*(1-s)+ s*((lam-Q)./(2*P));
        end
        
        function [W,ncflag]=calcVolRashid3dIter(C,c,bc,detC,v,vol,vi,i,Ielem2)
            ncflag=0;
            err=1e-5;
            
            %vols=zeros(size(C,1),1); %% until I figure out better
            %s=1;  % until I figure out better
   
            va=v(vi,:)-c(vi,:);
            
            P=(0.6*((.75/pi).^(2/3)))*(detC(vi).^(1/3));
            
            Q=sum(va.*VMath.SAx33(C(vi,:),va),2) + bc(vi);
            
            lam=(vol.*vol/3 +sum(Q./(2*P)))/sum(1./(2*P));
            %% need a newton loop here
            %lam=max(Q);
            f=99;
            fp=1;
            ii=1;
            
            while abs(f/fp/lam)>err&& ii<20,
                ker=(lam-Q)./(5*P/3);
                f=vol-sum(sqrt(ker).^3);
                fp=-1.5*sum(sqrt(ker)./(5*P/3));
                lam=lam-f/fp;
                if fp==0,
                    error('Derivative is equal to zero')
                end
                ii=ii+1;
            end
            if ii>=20,
                if abs(imag(f/fp/lam))>err,
                    if nargout<2,
                        warning('Too many Newton iterations')
                    else
                        ncflag=1;
                    end
                end
            end
            %%%
            
            W=sqrt((lam-Q)./(5*P/3)).^3;
            
            % [min(real(W)) max(real(W))]
            %vols=vols(vi)*(1-s)+ s*W;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%
        function [I1,I2]=setdiff_cell(I1,I2)
            for i=1:length(I2),
                vv=I2{i};
                for j=1:length(vv),
                    ind=find(I1{I2{i}(j)}==i,1);
                    if isempty(ind),
                        disp('hi')
                        vv(j)=[];
                    end
                end
                I2{i}=vv;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%
        function plot(ffrom,fto,fromelem,toelem)
            h=ffrom.plote(fromelem,[1 0 0]);
            fto.plote(toelem,[0 0 0],h);
            title(sprintf('%d From Elements(red), %d To Elements (black)',length(fromelem),length(toelem)))
        end
    end
end
