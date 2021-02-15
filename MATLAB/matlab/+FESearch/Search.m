classdef Search
    properties
        Tol = 0.01
        MaxNpts = 10e6
    end
    methods
        function obj=Search(Tol)
            if nargin>0,
                if ~isempty(Tol),
                    obj.Tol=Tol;
                end
            end
        end
        function [I]=ParticleSearch(obj,nodexyz,conn,findxyz,dxyz)
            %
            %
            %nodexyz -> master nodal coordinates
            %
            %conn -> a matrix of nodes of the master surface
            %
            %  findxyz -> slave nodes to be bounded
            %  dxyz -> 1x3 row vector maximum displacement during timestep of the
            %          slave nodes.  dxyz_max=[max(abs(dx)) max(abs(dy)) max(abs(dz))];
            %
            %  [I]=obj.ParticleSearch(nodexyz,conn,findxyz,dxyz)
            %
            %   Length of I is the number of elements in conn
            %   Each cell of I contains the node numbers in findxyz that potentially
            %   fall in the element I{i}
            %
            %%%%
            nc=size(conn,1);
            %nf=size(findxyz,1);
            %%%%
            xmin=min(reshape(nodexyz(conn,1),size(conn)),[],2)-(dxyz(1)+obj.Tol);
            xmax=max(reshape(nodexyz(conn,1),size(conn)),[],2)+(dxyz(1)+obj.Tol);
            ymin=min(reshape(nodexyz(conn,2),size(conn)),[],2)-(dxyz(2)+obj.Tol);
            ymax=max(reshape(nodexyz(conn,2),size(conn)),[],2)+(dxyz(2)+obj.Tol);
            zmin=min(reshape(nodexyz(conn,3),size(conn)),[],2)-(dxyz(3)+obj.Tol);
            zmax=max(reshape(nodexyz(conn,3),size(conn)),[],2)+(dxyz(3)+obj.Tol);

            I=cell(nc,1);
            %%%%%%
            %[Rx,Ix]=sort([xmin;xmax;findxyz(:,1)]);[ju,Ixp]=sort(Ix);Ix(Ix<2*nc)=
            %[Ry,Iy]=sort([ymin;ymax;findxyz(:,2)]);[ju,Iyp]=sort(Iy);
            %[Rz,Iz]=sort([zmin;zmax;findxyz(:,3)]);[ju,Izp]=sort(Iz);

            %for i=1:nc,
            %    ix=intersect(Ix(Ixp(i)+1:Ixp(i+nc)-1),intersect(Iy(Iyp(i)+1:Iyp(i+nc)-1),Iz(Izp(i)+1:Izp(i+nc)-1)));
            %    I{i}=ix(ix>2*nc)-2*nc;
            %end
            %%%%%%

            for i=1:nc,
                I{i}=find((xmin(i)<=findxyz(:,1)) & (xmax(i)>=findxyz(:,1)) & ...
                    (ymin(i)<=findxyz(:,2)) & (ymax(i)>=findxyz(:,2)) & ...
                    (zmin(i)<=findxyz(:,3)) & (zmax(i)>=findxyz(:,3)));
                %I{i}=setdiff(ix,conn(i,:));
            end
        end
        function [node,elem]=ElemID(obj,nodexyz,conn,findxyz,elemtype,dxyz)

            [I]=obj.ParticleSearch(nodexyz,conn,findxyz,dxyz);

            %% invert I
            %Iinv=obj.InvertSearch(I);
            [node,elem]=obj.InvertSearch2(I,elemtype,conn,nodexyz,findxyz);
        end
        function [I]=ElementSearch(obj,nodexyz,conn,findxyz,findconn,dxyz)
            %%% Identifies which elements in findxyz/findconn potentially are in
            %%%    nodexyz/conn

            nc=size(conn,1);
            %%% find the boundaries of the elements
            xmin=min(reshape(nodexyz(conn,1),size(conn)),[],2)-(dxyz(1)+obj.Tol);
            xmax=max(reshape(nodexyz(conn,1),size(conn)),[],2)+(dxyz(1)+obj.Tol);
            ymin=min(reshape(nodexyz(conn,2),size(conn)),[],2)-(dxyz(2)+obj.Tol);
            ymax=max(reshape(nodexyz(conn,2),size(conn)),[],2)+(dxyz(2)+obj.Tol);
            zmin=min(reshape(nodexyz(conn,3),size(conn)),[],2)-(dxyz(3)+obj.Tol);
            zmax=max(reshape(nodexyz(conn,3),size(conn)),[],2)+(dxyz(3)+obj.Tol);

            fxmin=min(reshape(findxyz(findconn,1),size(findconn)),[],2)-(dxyz(1)+obj.Tol);
            fxmax=max(reshape(findxyz(findconn,1),size(findconn)),[],2)+(dxyz(1)+obj.Tol);
            fymin=min(reshape(findxyz(findconn,2),size(findconn)),[],2)-(dxyz(2)+obj.Tol);
            fymax=max(reshape(findxyz(findconn,2),size(findconn)),[],2)+(dxyz(2)+obj.Tol);
            fzmin=min(reshape(findxyz(findconn,3),size(findconn)),[],2)-(dxyz(3)+obj.Tol);
            fzmax=max(reshape(findxyz(findconn,3),size(findconn)),[],2)+(dxyz(3)+obj.Tol);
            %%%
            I=cell(nc,1);
            for i=1:nc,
                I{i}=find(((xmin(i) >= fxmin) & (xmin(i) <= fxmax) | ...
                    (xmax(i) >= fxmin) & (xmax(i) <= fxmax) | ...
                    (xmin(i) <= fxmin) & (xmax(i) >= fxmax)) & ...
                    ((ymin(i) >= fymin) & (ymin(i) <= fymax) | ...
                    (ymax(i) >= fymin) & (ymax(i) <= fymax) | ...
                    (ymin(i) <= fymin) & (ymax(i) >= fymax)) & ...
                    ((zmin(i) >= fzmin) & (zmin(i) <= fzmax) | ...
                    (zmax(i) >= fzmin) & (zmax(i) <= fzmax) | ...
                    (zmin(i) <= fzmin) & (zmax(i) >= fzmax)));
            end

        end
        function [I]=BlockElementSearch(obj,nodexyz,blocks,findxyz,findblks,dxyz)
            %%% Identifies which elements in findxyz/findconn potentially are in
            %%%    nodexyz/conn, blocks is an exodus block format

            xmin=[];
            xmax=[];
            ymin=[];
            ymax=[];
            zmin=[];
            zmax=[];

            %%%
            for i=1:length(blocks),
                sz=size(blocks(i).Connectivity);
                xmin=cat(1,xmin,min(reshape(nodexyz(blocks(i).Connectivity,1),sz),[],2)-(dxyz(1)+obj.Tol));
                xmax=cat(1,xmax,max(reshape(nodexyz(blocks(i).Connectivity,1),size(blocks(i).Connectivity)),[],2)+(dxyz(1)+obj.Tol));
                ymin=cat(1,ymin,min(reshape(nodexyz(blocks(i).Connectivity,2),size(blocks(i).Connectivity)),[],2)-(dxyz(2)+obj.Tol));
                ymax=cat(1,ymax,max(reshape(nodexyz(blocks(i).Connectivity,2),size(blocks(i).Connectivity)),[],2)+(dxyz(2)+obj.Tol));
                zmin=cat(1,zmin,min(reshape(nodexyz(blocks(i).Connectivity,3),size(blocks(i).Connectivity)),[],2)-(dxyz(3)+obj.Tol));
                zmax=cat(1,zmax,max(reshape(nodexyz(blocks(i).Connectivity,3),size(blocks(i).Connectivity)),[],2)+(dxyz(3)+obj.Tol));
            end

            fxmin=[];
            fxmax=[];
            fymin=[];
            fymax=[];
            fzmin=[];
            fzmax=[];

            for i=1:length(findblks),
                sz=size(findblks(i).Connectivity);
                fxmin=cat(1,fxmin,min(reshape(findxyz(findblks(i).Connectivity,1),sz),[],2)-(dxyz(1)+obj.Tol));
                fxmax=cat(1,fxmax,max(reshape(findxyz(findblks(i).Connectivity,1),sz),[],2)+(dxyz(1)+obj.Tol));
                fymin=cat(1,fymin,min(reshape(findxyz(findblks(i).Connectivity,2),sz),[],2)-(dxyz(2)+obj.Tol));
                fymax=cat(1,fymax,max(reshape(findxyz(findblks(i).Connectivity,2),sz),[],2)+(dxyz(2)+obj.Tol));
                fzmin=cat(1,fzmin,min(reshape(findxyz(findblks(i).Connectivity,3),sz),[],2)-(dxyz(3)+obj.Tol));
                fzmax=cat(1,fzmax,max(reshape(findxyz(findblks(i).Connectivity,3),sz),[],2)+(dxyz(3)+obj.Tol));
            end
            %%%
            nc=length(xmin);
            I=cell(nc,1);
            for i=1:nc,
                I{i}=find(((xmin(i) >= fxmin) & (xmin(i) <= fxmax) | ...
                    (xmax(i) >= fxmin) & (xmax(i) <= fxmax) | ...
                    (xmin(i) <= fxmin) & (xmax(i) >= fxmax)) & ...
                    ((ymin(i) >= fymin) & (ymin(i) <= fymax) | ...
                    (ymax(i) >= fymin) & (ymax(i) <= fymax) | ...
                    (ymin(i) <= fymin) & (ymax(i) >= fymax)) & ...
                    ((zmin(i) >= fzmin) & (zmin(i) <= fzmax) | ...
                    (zmax(i) >= fzmin) & (zmax(i) <= fzmax) | ...
                    (zmin(i) <= fzmin) & (zmax(i) >= fzmax)));
            end
        end
        function [I]=BlockSearch(obj,nodexyz,blocks,findxyz,findblks,dxyz)
            %%% Identifies which blocks in findxyz/findconn potentially
            %%%  overlapped with blocks.  blocks is an exodus block format

            xmin=zeros(length(blocks),1);
            xmax=zeros(length(blocks),1);
            ymin=zeros(length(blocks),1);
            ymax=zeros(length(blocks),1);
            zmin=zeros(length(blocks),1);
            zmax=zeros(length(blocks),1);

            %%%
            for i=1:length(blocks),
                sz=size(blocks(i).Connectivity);
                xmin(i)=min(min(reshape(nodexyz(blocks(i).Connectivity,1),sz),[],2)-(dxyz(1)+obj.Tol));
                xmax(i)=max(max(reshape(nodexyz(blocks(i).Connectivity,1),sz),[],2)+(dxyz(1)+obj.Tol));
                ymin(i)=min(min(reshape(nodexyz(blocks(i).Connectivity,2),sz),[],2)-(dxyz(2)+obj.Tol));
                ymax(i)=max(max(reshape(nodexyz(blocks(i).Connectivity,2),sz),[],2)+(dxyz(2)+obj.Tol));
                zmin(i)=min(min(reshape(nodexyz(blocks(i).Connectivity,3),sz),[],2)-(dxyz(3)+obj.Tol));
                zmax(i)=max(max(reshape(nodexyz(blocks(i).Connectivity,3),sz),[],2)+(dxyz(3)+obj.Tol));
            end

            fxmin=zeros(length(findblks),1);
            fxmax=zeros(length(findblks),1);
            fymin=zeros(length(findblks),1);
            fymax=zeros(length(findblks),1);
            fzmin=zeros(length(findblks),1);
            fzmax=zeros(length(findblks),1);

            for i=1:length(findblks),
                sz=size(findblks(i).Connectivity);
                fxmin(i)=min(min(reshape(findxyz(findblks(i).Connectivity,1),sz),[],2)-(dxyz(1)+obj.Tol));
                fxmax(i)=max(max(reshape(findxyz(findblks(i).Connectivity,1),sz),[],2)+(dxyz(1)+obj.Tol));
                fymin(i)=min(min(reshape(findxyz(findblks(i).Connectivity,2),sz),[],2)-(dxyz(2)+obj.Tol));
                fymax(i)=max(max(reshape(findxyz(findblks(i).Connectivity,2),sz),[],2)+(dxyz(2)+obj.Tol));
                fzmin(i)=min(min(reshape(findxyz(findblks(i).Connectivity,3),sz),[],2)-(dxyz(3)+obj.Tol));
                fzmax(i)=max(max(reshape(findxyz(findblks(i).Connectivity,3),sz),[],2)+(dxyz(3)+obj.Tol));
            end
            %%%
            nc=length(blocks);
            I=cell(nc,1);
            for i=1:nc,
                I{i}=find(((xmin(i) >= fxmin) & (xmin(i) <= fxmax) | ...
                    (xmax(i) >= fxmin) & (xmax(i) <= fxmax) | ...
                    (xmin(i) <= fxmin) & (xmax(i) >= fxmax)) & ...
                    ((ymin(i) >= fymin) & (ymin(i) <= fymax) | ...
                    (ymax(i) >= fymin) & (ymax(i) <= fymax) | ...
                    (ymin(i) <= fymin) & (ymax(i) >= fymax)) & ...
                    ((zmin(i) >= fzmin) & (zmin(i) <= fzmax) | ...
                    (zmax(i) >= fzmin) & (zmax(i) <= fzmax) | ...
                    (zmin(i) <= fzmin) & (zmax(i) >= fzmax)));
            end
        end
        function [I]=NearestPoint(obj,nodexyz,findxyz,type)
            %
            % this function finds the nearest point in list findxyz
            %   to each point in nodexyz and returns the indices
            %  I is the length of findxyz, ie each point of findxyz will
            %  have a index from nodexyz which is closest to the point in
            %  findxyz.
            %
            if nargin<4,
                if size(nodexyz,1)>1e3 && size(findxyz,1)>1e3,
                    type=2;
                else
                    type=1;
                end
            end
            %%%
            switch type
                case 1
                    I=dsearchn(nodexyz,findxyz);
                case 2
                    I=zeros(size(findxyz,1),1);
                    reterr=zeros(size(findxyz,1),1);

                    for i=1:size(findxyz,1),
                        err=(nodexyz(:,1)-findxyz(i,1)).^2 + (nodexyz(:,2)-findxyz(i,2)).^2;
                        if size(nodexyz,2)==3,
                            err=err+(nodexyz(:,3)-findxyz(i,3)).^2;
                        end
                        [reterr(i),I(i)]=min(err);
                    end
                    % [idx,I]=sort(idx);
                    % reterr=reterr(I);

                otherwise
                    error('Unknown search option')
            end
        end
        %
    end
    methods (Static)
        function Iinv=InvertSearch(I)
            % figure out how many elements you are dealing with
            n=0;
            for i=1:length(I),
                n=n+length(I{i});
            end
            % Allocate a new I
            nodes=zeros(n,1);
            elem=zeros(n,1);
            % need to redefine indexing for for nodes(n*(i-1):n*i) or something
            for i=1:length(I),
                if ~isempty(I{i}),
                    n=length(I{i});
                    nodes((n*(i-1)+1):(n*i))=I{i};
                    elem((n*(i-1)+1):(n*i))=i*ones(n,1);
                end
            end
            un=unique(nodes);
            n=length(un);
            Iinv=cell(n,1);
            for i=1:n,
                Iinv{i}=elem(un(i)==nodes);
            end

        end
        function [nodes,elems]=InvertSearch2(I,elemtype,conn,nodexyz,findxyz)
            % figure out how many elements you are dealing with
            n=0;
            for i=1:length(I),
                if ~isempty(I{i})
                    n=n+length(I{i});
                end
            end
            % Allocate a new I
            node=zeros(n,1);
            elem=zeros(n,1);
            ii=1;
            % need to redefine indexing for for nodes(n*(i-1):n*i) or something
            for i=1:length(I),
                if ~isempty(I{i}),
                    n=length(I{i});
                    if (any(~I{i})),
                        disp(num2str(i))
                    end
                    node(ii:(ii+n-1))=I{i};
                    elem(ii:(ii+n-1))=i*ones(n,1);
                    ii=ii+n;
                end
            end

            un=unique(node);
            n=length(un);
            %Iinv=cell(n,1);
            nodes=1:n;
            elems=zeros(n,1);
            Sobj=ShapeFactory.CreateShape(elemtype,conn,nodexyz);
            for i=1:n,
                iii=find(un(i)==node);

                Iinv=elem(un(i)==node);

                if ~isempty(Iinv),
                    pts=Sobj.Global2Local(repmat(findxyz(i,:),length(Iinv),1),Iinv);
                    if any(isnan(pts)),
                        findxyz_nan=findxyz(i,:);
                        Iinv_nan=Iinv;
                        error('NANS were found in InvertSearch2')
                    end

                    [ju,idx]=min(sum(pts.*pts,2));

                    elems(i)=Iinv(idx);
                end
            end

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
            %max(localxyz(idx2,:),[],1)
            %min(localxyz(idx2,:),[],1)
        end
    end
end
