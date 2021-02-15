classdef ElementNode
    properties
        ModelObj
        BlocksOfInterest
        Tol = 0.6  % Amount to allow for extrapolation
        SearchTol = .05  % use default in the search
        TightTol = 0.5 % amount where it doesn't matter whether you extrapolate or interpolate
        MaxLength = 1e6 % maximum number of local points to evaluate at a time to save memory.
    end
    methods
        function obj=ElementNode(fexo,idxblk,tol,stol)
            if nargin>0,
                obj.ModelObj=fexo;
                obj.BlocksOfInterest=idxblk;
                if nargin>2 && ~isempty(tol),
                    obj.Tol=tol;
                end
                if nargin>3 && ~isempty(stol),
                    obj.SearchTol=stol;
                end
            end
        end
        function [lcoords,elems,blks,nodes,nodesnotfound]=LocalCoord(obj,locatenodes)
            %% elem(:), node(:) and globalxyz all have to have the same number of
            %% rows.
            
            searchObj=FESearch.Search(obj.SearchTol);
            elems=zeros(size(locatenodes,1),1);
            blks=zeros(size(locatenodes,1),1);
            nodes=zeros(size(locatenodes,1),1);
            nodesleft= true(size(locatenodes,1),1);
            locatenodenum=1:size(locatenodes,1);
            lcoords=zeros(size(locatenodes,1),3);
            jj=1;
            for i=1:length(obj.BlocksOfInterest),
                blkid=obj.BlocksOfInterest(i);
                [nodexyz,conn]=obj.global2local_conn(obj.ModelObj.Nodes.Coordinates, ...
                    obj.ModelObj.Blocks(blkid).Connectivity);
                %%%%%%%%%%%%%%%Remove
                if 1,
                    I=searchObj.ParticleSearch(nodexyz,conn,locatenodes(nodesleft,:),[0 0 0]);
                    %save data I
                else
                    %load data
                end
                nnodes=0;
                for j=1:length(I),
                    %%This loop figures out how long to make elem
                    if ~isempty(I{j}),
                        nnodes=nnodes+length(I{j});
                    end
                end
                
                elem=zeros(nnodes,1);
                node=zeros(nnodes,1);
                blk=zeros(nnodes,1);
                ii=1;
                nodenum=locatenodenum(nodesleft);
                for j=1:length(I),
                    if ~isempty(I{j}),
                        elem(ii:(ii+length(I{j})-1))=j;
                        node(ii:(ii+length(I{j})-1))=nodenum(I{j}(:));
                        blk(ii:(ii+length(I{j})-1))=i;
                        ii=ii+length(I{j});
                    end
                end
                %% Calculate the local coordinates
                [nodexyz,conn]=FESearch.ElementNode.global2local_conn(obj.ModelObj.Nodes.Coordinates,obj.ModelObj.Blocks(blkid).Connectivity);
                shapeObj=ShapeFactory.CreateShape(obj.ModelObj.Blocks(blkid).ElementType, ...
                    conn, nodexyz);
               
                maxlen=obj.MaxLength;
                nnn=length(node);
                ind=1:nnn;
                
                nl=ceil(nnn/maxlen);
                pts=zeros(0,3);
                nodehold=zeros(0,1);
                elemhold=zeros(0,1);
                blkhold=zeros(0,1);
                for k=1:nl,
                    
                    if k<nl,
                        vii=ind(((k-1)*maxlen+1):(k*maxlen));
                    else
                        vii=ind(((k-1)*maxlen+1):end);
                    end
                    ptsh=shapeObj.Global2Local(locatenodes(node(vii),:),elem(vii));
                    if size(pts,2)<3,
                        ptsh=[ptsh zeros(size(pts,1),(3-size(pts,2)))];
                    end
                    if any(isnan(ptsh)),
                        fprintf('There are %d NAN entries\n',length(find(isnan(ptsh(:,1)))))
                    end
                    idx=FESearch.ElementNode.IsValid(shapeObj,ptsh,node(vii),obj.Tol,obj.TightTol);
                    
                    nodehold=cat(1,nodehold,node(vii(idx)));
                    elemhold=cat(1,elemhold,elem(vii(idx)));
                    blkhold=cat(1,blkhold,blk(vii(idx)));
                    pts=cat(1,pts,ptsh(idx,:));
                end
                
                idx=FESearch.ElementNode.IsValid(shapeObj,pts,nodehold,obj.Tol,obj.TightTol);
                
                elems(jj:jj+length(idx)-1)=elemhold(idx);
                blks(jj:jj+length(idx)-1)=blkhold(idx);
                nodes(jj:jj+length(idx)-1)=nodehold(idx);
                lcoords(jj:jj+length(idx)-1,:)=pts(idx,:);
                jj=jj+length(idx);
                nodesleft(nodehold(idx))=false;
            end
            a=1:size(locatenodes,1);
            nodesnotfound=a(nodesleft);
            [ju,I]=sort(nodes);
            idx=length(nodesnotfound)+1:length(nodes);
            %%
            nodes=nodes(I(idx));
            elems=elems(I(idx));
            blks=blks(I(idx));
            lcoords=lcoords(I(idx),:);
        end
    end
    methods (Static)
        function [I,idx]=NonRepeated(In)
            [Insort,idx1]=sort(In);
            [ju,idx]=setdiff(In,In(idx1(~diff(Insort))));
            idx=sort(idx);
            I=In(idx);
        end
        function [I]=Repeated(In)
            [Insort]=sort(In);
            I=Insort(~diff(Insort));
        end
        function idx2=IsValid(ShapeObj,localxyz,node,tol,tighttol)
            [dx,dy,dz]=ShapeObj.LocalBounds;
            if size(localxyz,2)>2,
                idx=find(~isnan(localxyz(:,1)) & ...
                    (localxyz(:,1)>=(dx(1)-tol) & localxyz(:,1)<=(dx(2)+tol)) ...
                    &(localxyz(:,2)>=(dy(1)-tol) & localxyz(:,2)<=(dy(2)+tol)) ...
                    & (localxyz(:,3)>=(dz(1)-tol) & localxyz(:,3)<=(dz(2)+tol)));
            elseif size(localxyz,2)>1,
                idx=find(~isnan(localxyz(:,1)) & ...
                    (localxyz(:,1)>=(dx(1)-tol) & localxyz(:,1)<=(dx(2)+tol)) ...
                    &(localxyz(:,2)>=(dy(1)-tol) & localxyz(:,2)<=(dy(2)+tol)));
            else
                idx=find(~isnan(localxyz(:,1)) & ...
                    (localxyz(:,1)>=(dx(1)-tol) & localxyz(:,1)<=(dx(2)+tol)));
            end

            switch 2
                case 1 % this just picks the first node as the good one.  Not optimal but fast
                    [node,ii]=unique(node(idx));
                    idx2=idx(ii);
                case 2
                    nodechk=FESearch.ElementNode.Repeated(node(idx));
                    if ~isempty(nodechk),
                        if size(localxyz,2)>2,
                            idx2=find((localxyz(:,1)>=(dx(1)-tighttol) & localxyz(:,1)<=(dx(2)+tighttol)) ...
                                &(localxyz(:,2)>=(dy(1)-tighttol) & localxyz(:,2)<=(dy(2)+tighttol)) ...
                                & (localxyz(:,3)>=(dz(1)-tighttol) & localxyz(:,3)<=(dz(2)+tighttol)));
                        else
                            idx2=find((localxyz(:,1)>=(dx(1)-tighttol) & localxyz(:,1)<=(dx(2)+tighttol)) ...
                                &(localxyz(:,2)>=(dy(1)-tighttol) & localxyz(:,2)<=(dy(2)+tighttol)));
                        end
                        [node2,ii]=unique(node(idx2));
                        idx2=idx2(ii);

                        %
                        [node3,idx3]=setdiff(node(idx),node(idx2));
                        % if they are within tighttol
                        %%  then it doesn't matter which element they are in
                        %%%%%%%
                        %% This code looks at the two values and chooses the one
                        %% with the smallest maximum local coordinate.  Probably
                        %% optimal but very slow.  This code needs to be optimized.
                        for i=1:length(node3),
                            ind=find(node3(i)==node(idx));
                            if length(ind)>1,
                                [ju]=max(abs(localxyz(idx(ind),:)),[],2);
                                [ju1,idxhold]=min(ju);
                                idx2=cat(1,idx2,idx(ind(idxhold)));
                            else
                                idx2=cat(1,idx2,idx(ind));
                            end
                        end
                    else
                        idx2=idx;
                    end
                case 3  % an optimal version of case 2
                    % for locatenodes that fall within multiple elements,
                    % use the local coordinates that have the smallest
                    % absolute max value.
                    %[localxyz(idx,:) node(idx)]
                    [lmax]=max(abs(localxyz(idx,:)),[],2);
                    [lmax,visort]=sort(lmax);
                    %[lmax node(idx(visort))]
                    [node,ii]=unique(node(idx(visort)),'first');
                    %[lmax(ii) node]
                    idx2=idx(visort(ii));
                    %localxyz(idx2,:)
                    %length(node(idx))
            end
            if 0,
                max(localxyz(idx2,:),[],1)
                min(localxyz(idx2,:),[],1)
            end
            
        end
        function [nodexyz_local,conn_local]=global2local_conn(nodexyz,conn)
            
            idx=unique(conn(:));
            nodexyz_local=nodexyz(idx,:);
            %
            a=zeros(size(nodexyz_local,1),1);
            a(idx)=(1:length(idx))';
            conn_local=a(conn);
        end
    end
end
