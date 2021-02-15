classdef ElemInterp < FEMap.Map
    properties
        AbsoluteSearchTolerance = []
        LocalCoordinateTolerance = []
    end
    methods
        function obj=ElemInterp(fexoFrom,fexoTo)
            if nargin>0,
                obj.ModelObjTo=fexoTo;
                obj.ModelObjFrom=fexoFrom; 
            end
        end
        function obj=E2EMap(obj,blkidxfrom,blkidxto,numgausspts,enamefrom,enameto,tidx,type)
            ngauss=1;
            switch type
                case 1  % use a virtual mesh connecting centroids to map
                    obj=VirtualMeshInterpElem(obj,enamefrom,numgausspts,tidx);
                case 2 % find the closest centroids to map variables
                    [ptxyzFrom,evarFrom]=obj.getEVar(blkidxfrom,enamefrom,tidx);
                    [enodexyz,nelem,gp_id]=obj.CalcElemCoords(obj.ModelObjTo, ...
                        blkidxto,numgausspts);
                    Sobj=FESearch.Search;
                    elem=Sobj.NearestPoint(ptxyzFrom,enodexyz,2);
                    nodes=1:size(enodexyz,1);
  
                    evar=zeros(size(enodexyz,1),1);
                    evar(nodes)=evarFrom(elem);
                    obj=obj.WriteElem(evar,enameto, ...
                        blkidxto,nelem,ngauss,gp_id,obj.ModelObjFrom.Time(tidx));
                case 3 % determine if an element centroid lies within an element then directly map
                    [enodexyz,nelem,gp_id]=obj.CalcElemCoords(obj.ModelObjTo, ...
                        blkidxto,numgausspts);
                    %%%
                    [ptxyzFrom,evarFrom]=obj.getEVar(blkidxfrom,enamefrom,tidx);
                    %%%%
                    SrchObj=FESearch.ElementNode(obj.ModelObjFrom,blkidxfrom);
                    %%%
                    [lc,elem,blks, ...
                        points,pointsnotmapped]=SrchObj.LocalCoord(enodexyz);
                    if ~isempty(pointsnotmapped),
                        warning('%d points not mapped',length(pointsnotmapped));
                    end
                    
                    nodes=points;
                    evar=zeros(size(enodexyz,1),1);
                    evar(nodes)=evarFrom(elem);
                    obj=obj.WriteElem(evar,enameto, ...
                        blkidxto,nelem,ngauss,gp_id,obj.ModelObjFrom.Time(tidx));
                case 4  % development. should do the same as 3...
                    [enodexyz,nelem,gp_id]=obj.CalcElemCoords(obj.ModelObjTo, ...
                        blkidxto,numgausspts);
                    
                    %%%
                    [ptxyzFrom,evarFrom]=obj.getEVar(blkidxfrom,enamefrom,tidx);
                    %%%%
                    Sobj=FESearch.Search;
                    nodes=[];
                    elem=[];
                    for i=1:length(blkidxfrom),
                        [no,el]=Sobj.ElemID(obj.ModelObjFrom.Nodes.Coordinates, ...
                            obj.ModelObjFrom.Blocks(blkidxfrom(i)).Connectivity,enodexyz, ...
                            obj.ModelObjFrom.Blocks(blkidxfrom(i)).ElementType,[0 0 0]);
                        nodes=[nodes;no(:)];
                        elem=[elem;el(:)];
                    end
                    evar=zeros(size(enodexyz,1),1);
                    evar(nodes)=evarFrom(elem);
                    obj=obj.WriteElem(evar,enameto, ...
                        blkidxto,nelem,ngauss,gp_id,obj.ModelObjFrom.Time(tidx));
                case 5 % just calculate the total energy and put it uniformly on the To mesh
                    tot_from_evar=zeros(length(blkidxfrom),1);
                    tot_vol=zeros(length(blkidxfrom),1);
                    for j=1:length(blkidxfrom), % figure out total energy on original mesh
                        [tot_from_evar(j),tot_vol(j)]=FEMap.Map.IntBlock(obj.ModelObjFrom,blkidxfrom(j), ...
                            enamefrom,tidx);
                    end
                    [enodexyz,nelem,gp_id]=obj.CalcElemCoords(obj.ModelObjTo, ...
                        blkidxto,numgausspts);
                    evar=zeros(size(enodexyz,1),1);
                    vol=obj.ModelObjTo.Volume(blkidxto);
                    for i=1:length(blkidxto),
                        evar(1:(nelem(i)*size(gp_id,1)))=sum(tot_from_evar)/vol;
                    end
                    obj=obj.WriteElem(evar,enameto, ...
                        blkidxto,nelem,ngauss,gp_id,obj.ModelObjFrom.Time(tidx));
            end
            
            obj.VariableTypeTo='Element';
            obj.VariableTypeFrom='Element';
        end
        function obj=E2NMap(obj,blkidxfrom,enamefrom,nnameto,type,tidx)
            
            switch type
                case 1 %% Nodal Averaging
                    names=obj.ModelObjFrom.ElemVars.getNames;
                    evaridx=find(strncmp(enamefrom,names,length(enamefrom)));
                    if isempty(evaridx),
                        error('Element Variable name not found')
                    end
                    [nodal_data]=obj.NodalAverage(blkidxfrom,evaridx,tidx);
                case 2 %% Local Least Squares
                    names=obj.ModelObjFrom.ElemVars.getNames;
                    %evaridx=strmatch(ename,names);
                    evaridx=find(strncmp(enamefrom,names,length(enamefrom)));
                    if isempty(evaridx),
                        error('Element Variable name not found')
                    end
                    [nodal_data]=obj.LocalLeastSquares(blkidxfrom,evaridx,tidx);
                case 3
                    [ptxyzFrom,evarFrom]=obj.getEVar(blkidxfrom,enamefrom,tidx);
                    nodal_data=griddatan(ptxyzFrom,evarFrom,obj.ModelObjTo.Nodes.Coordinates,'linear',{'Qt','Qbb','Qc','Qz'});
                    if any(isnan(nodal_data)),
                        error('Mesh of element coordinates does not encompass the nodal points. Try type=4.');
                    end
                    %%  TODO: Need to address NaN's (points not internal to
                    %%  ptxyz)  actually this will always fail since the
                    %%  nodes are outside the elements duh.
                case 4 % this is the same as 3 except that it uses my code for the interpolation which allows for extrapolation
                    nodal_data=obj.VirtualMeshInterpNodal(enamefrom,tidx); 
            end
            time=obj.ModelObjFrom.Time(tidx);
            
            obj.ModelObjFrom=obj.ModelObjFrom.AddNodalVar(nnameto,nodal_data,time);
            obj.VariableTypeTo='Nodal';
            obj.VariableTypeFrom='Element';
        end
    end
    methods (Access=private)
        function [nodal_data]=NodalAverage(obj,blkidx,evaridx,tidx)
            % data is a elem x 1
            nodal_data=zeros(size(obj.ModelObjFrom.Nodes.Coordinates,1),1);
            node_used=zeros(size(nodal_data));
            for i=1:length(blkidx),
                [node_elem,nodes]=inv_connect(obj,blkidx);
                elem_data=obj.ModelObjFrom.ElemVars(blkidx(i),evaridx).Data(:,tidx);
                for j=1:length(nodes),
                    nodal_data(nodes(j))=nodal_data(nodes(j))+mean(elem_data(node_elem{j}));
                    node_used(nodes(j))=node_used(nodes(j))+1;
                end
            end
            if any(node_used>1),
                nodal_data(node_used)=nodal_data(node_used)./node_used(node_used);
            end
        end
        function [nodal_data]=LocalLeastSquares(obj,blkidx,evaridx,tidx)
            % data is a elem x 1
            nodal_data=zeros(size(obj.ModelObjFrom.Nodes.Coordinates,1),1);
            node_used=zeros(size(nodal_data));
            mean_nodes=zeros(size(nodal_data));
            for i=1:length(blkidx),
                [node_elem,nodes]=inv_connect(obj,blkidx);
                elem_cent=FEMap.ElemInterp.elemcentroid(obj.ModelObjFrom.Nodes.Coordinates,obj.ModelObjFrom.Blocks(blkidx).Connectivity);
                elem_data=obj.ModelObjFrom.ElemVars(blkidx(i),evaridx).Data(:,tidx);

                for j=1:length(nodes),
                    [x,r]=Elem2Nodal.pinv([ones(length(node_elem{j}),1) elem_cent(node_elem{j},:)], ...
                        elem_data(node_elem{j}));
                    if r>3,
                        nodal_data(nodes(j))=nodal_data(nodes(j))+[1 obj.ModelObjFrom.Nodes.Coordinates(nodes(j),:)]*x;
                    else
                        nodal_data(nodes(j))=mean(elem_data(node_elem{j}));
                        mean_nodes(nodes(j))=1;            
                    end
                    node_used(nodes(j))=node_used(nodes(j))+1;
                end
                % Heuristics for boundary nodes
                idx=find(mean_nodes);
                for j=idx',
                    [nodal_data(j),nidx]=obj.bc1(blkidx(i),nodal_data,mean_nodes,node_elem{j});
                    %[nodal_data(j),nidx]=obj.bc2(blkidx(i),nodal_data,mean_nodes,node_elem{j},obj.ModelObjFrom.Nodes.Coordinates(j,:));
                    node_used(j)=node_used(nidx);
                end
                
            end
            if any(node_used>1),
                nodal_data(node_used)=nodal_data(node_used)./node_used(node_used);
            end
        end
        
        function [nodal_data,nodes,node_elem]=LocalLeastSquares_old(obj,elem_data)
            [node_elem,nodes]=inv_connect(obj,blkidx);
            elem_cent=Elem2Nodal.elemcentroid(obj.ModelObjFrom.Nodes.Coordinates,obj.Blocks(blkidx).Connectivity);
            nodal_data=zeros(size(elem_data,1),length(nodes));
            mean_nodes=zeros(length(nodes),1);
            for i=1:length(nodes),
                [x,r]=Elem2Nodal.pinv([ones(length(node_elem{i}),1) elem_cent(node_elem{i},:)], ...
                    elem_data(:,node_elem{i})');
                if r>3,
                    nodal_data(:,i)=[1 obj.ModelObjFrom.Nodes.Coordinates(nodes(i),:)]*x;
                else
                    nodal_data(:,i)=mean(elem_data(:,node_elem{i}),2)';
                    mean_nodes(i)=1;
                end
            end
            % Heuristics for boundary nodes 
            idx=find(mean_nodes);
            for i=idx',
                %nodal_data(:,i)=obj.bc1(nodal_data,mean_nodes,node_elem{i});
                nodal_data(:,i)=obj.bc2(nodal_data,mean_nodes,node_elem{i},obj.ModelObjFrom.Nodes.Coordinates(i,:));
            end
                
        end
     
        function [node_elem,nodes]=inv_connect(obj,blkidx)
            %
            %
            num_nodes_per_element=size(obj.ModelObjFrom.Blocks(blkidx).Connectivity,2);
            %
            conn=obj.ModelObjFrom.Blocks(blkidx).Connectivity';
            [scon,idx]=sort(conn(:));

            [nodes,idx2]=unique(scon,'last');
            idx=idx-1;
            
            idx=floor(idx/num_nodes_per_element)+1; % these are the elements corresponding to the nodes in scon
            idx(~idx)=1;
            %
            node_elem=cell(length(nodes),1);
            idxprev=1;
            for i=1:length(nodes),
                node_elem{i}=idx(idxprev:idx2(i));
                idxprev=idx2(i)+1;
            end
        end
        function [nd,node]=bc1(obj,blkidx,nodal_data,mean_nodes,elems)
            % I think this uses the maximum value in the same elements as
            % the node of interest
            enode=obj.ModelObjFrom.Blocks(blkidx).Connectivity(elems,:);
            nodes=unique(enode(:));
            local=nodal_data(nodes);
            local=local(~mean_nodes(nodes));
            nodes=nodes(~mean_nodes(nodes));
            [nd,i]=max(abs(local));
            node=nodes(i);
        end
        function [nd,node]=bc2(obj,blkidx,nodal_data,mean_nodes,elems,nodexyz)
            %% I think this uses the nearest value to the point of interest
            %% nodexyz is the coord of node of interest
            enode=obj.ModelObjFrom.Blocks(blkidx).Connectivity(elems,:);
            nodes=unique(enode(:));
            local=nodal_data(nodes);
            local=local(~mean_nodes(nodes));
            coord=obj.ModelObjFrom.Nodes.Coordinates(nodes,:);
            coord=coord(~mean_nodes(nodes),:);
            nodes=nodes(~mean_nodes(nodes));
            dist=sum((coord-repmat(nodexyz,size(coord,1),1)).*(coord-repmat(nodexyz,size(coord,1),1)),2);
            [ju,i]=min(dist);
            nd=local(i);
            node=nodes(i);
        end
        function [nodal_data]=VirtualMeshInterpNodal(obj,ename,tidx)
            [ptxyzFrom,evarFrom]=obj.getEVar(blkidxfrom,ename,tidx);
            T=delaunayn(ptxyzFrom,{'Qt','Qbb','Qc','Qz'});
            fvirt=FEMesh.Exodus('tetra4',T,ptxyzFrom);
            SrchObj=FESearch.ElementNode(fvirt,1,obj.LocalCoordinateTolerance,obj.AbsoluteSearchTolerance);
            [lcoord,elem,blk, ...
                pts,ptsnotmapped]=SrchObj.LocalCoord(obj.ModelObjTo.Nodes.Coordinates);
            if ~isempty(ptsnotmapped),
                error('Need to increase search tolerances')
            end   
            
            shpobj=ShapeFactory.CreateShape('tet',T,ptxyzFrom);
            N=shpobj.Interpolate(lcoord);
            nodal_data=sum(N.*evarFrom(T(elem,:)),2);
        end
        function [obj]=VirtualMeshInterpElem(obj,ename,numgausspts,numintpts,tidx)
            [ptxyzFrom,evarFrom]=obj.getEVar(blkidxfrom,ename,tidx);
            T=delaunayn(ptxyzFrom,{'Qt','Qbb','Qc','Qz'});
            
            fvirt=FEMesh.Exodus('tetra4',T,ptxyzFrom);
            %%%
            [enodexyz,ngauss,nelem,gp_id,nint,weights,detJ]=obj.CalcElemCoords(obj.ModelObjTo, ...
                blkidxto,numgausspts,numintpts);
            
            %%%
            [elem,lcoord1]=tsearchn(ptxyzFrom,T,enodexyz);
            lcoord=lcoord1(:,1:3);
            idx=find(isnan(lcoord(:,1)));
            %%%
            if 0,
                fem=exo2imat3(fvirt);
                plotfem(fem)
                hold on
                plot3(enodexyz(:,1),enodexyz(:,2),enodexyz(:,3),'ko')
            end
            
            if ~isempty(idx),
                %% Extrapolation is done by just matching the points with
                %% the centroids of the elements.  Not perfect but fast
                SrchObj=FESearch.Search;
                ShpObj=ShapeFactory.CreateShape('tet',fvirt.Blocks(1).Connectivity,fvirt.Nodes.Coordinates);
                centxyz=ShpObj.GlobalCentroid;
                elem(idx)=SrchObj.NearestPoint(centxyz,enodexyz(idx,:));
                lcoord(idx,:)=ShpObj.Global2Local(enodexyz(idx,:),elem(idx));
%                 
            end
            %%%%
            shpobj=ShapeFactory.CreateShape('tet',T,ptxyzFrom);
            N=shpobj.Interpolate(lcoord);
            %%
            if 0,
                nodexyz=shpobj.Local2Global(lcoord,elem);
                plot3(nodexyz(:,1),nodexyz(:,2),nodexyz(:,3),'r.',enodexyz(:,1),enodexyz(:,2),enodexyz(:,3),'bo')
            end
            %%
            evar=sum(N.*evarFrom(T(elem,:)),2);
            obj=WriteElem(obj,evar,ename, ...
                        blkidxto,nelem,ngauss,gp_id,obj.ModelObjFrom.Time(tidx));
        end
        
    end
    methods (Static)
        function [x,r]=pinv(A,b)
            if size(A,1)<size(A,2),
                r=size(A,1);
                x=zeros(size(b));
                return
            end
            [u,ep,v]=svd(A);
            s=diag(ep);
            r=sum((s/s(1))>1e-8);
            x=v(:,1:r)*diag(ones(r,1)./s(1:r))*u(:,1:r)'*b;
        end
        function cent=elemcentroid(nodexyz,conn)
            %% this isn't correct
            n=size(nodexyz,1);    
            cent=[mean(nodexyz(conn),2) mean(nodexyz(n+conn),2) mean(nodexyz(2*n+conn),2)];
        end
        
    end
end