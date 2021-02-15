classdef Exodus < handle
    properties
        Title
        Filename
        FloatingPointWordSize  = 8
        FileSize               = 1
        LastTimeWritten        =[]
        QARecords
        InfoRecords
        Nodes
        CoordinateFrames
        Blocks
        ElementMap  =[]
        Nodesets
        Sidesets
        Nemesis
        SuperElement = []
        Time
        GlobalVars
        NodalVars
        ElemVars
        NodesetVars
        SidesetVars
    end
    methods
        function obj=Exodus(fpws,fs,coord)
            if nargin>0,
                
                if length(fpws)==1,
                    obj.FloatingPointWordSize=fpws;
                    obj.FileSize=fs;
                else
                    %% obj=Exodus(elemtype,conn,coord)
                    obj.Title='Default Matlab Title';
                    obj.Filename='junkppp.g';
                    a=ver('Matlab');
                    obj.QARecords=FEMesh.QARecords('Matlab',sprintf('Version: %s',a.Version));
                    obj.Nodes=FEMesh.Nodes(coord,[]);
                    obj.Blocks=FEMesh.Blocks(1);
                    obj.Blocks.ElementType=fpws;
                    obj.Blocks.Connectivity=fs;
                    obj.ElementMap=(1:size(fs,1))';
                    obj.Nodesets=FEMesh.Nodesets([]);
                    obj.Sidesets=FEMesh.Sidesets([]);
                    obj.GlobalVars=FEMesh.GlobalVars([]);
                    obj.NodalVars=FEMesh.NodalVars([]);
                    obj.ElemVars=FEMesh.ElemVars(1,[]);
                    obj.NodesetVars=FEMesh.NodesetVars([]);
                    obj.SidesetVars=FEMesh.SidesetVars([]);
                end
            end
        end
        function display(obj)
            fprintf('Database:  %s\n\n',obj.Filename)
            fprintf('%s\n\n',obj.Title)
            fprintf('Number of coordinates per node       = %13d\n',size(obj.Nodes.Coordinates,2))
            fprintf('Number of nodes                      = %13d\n',size(obj.Nodes.Coordinates,1))
            fprintf('Number of elements                   = %13d\n',obj.Blocks.numelem)
            fprintf('Number of element blocks             = %13d\n\n',length(obj.Blocks))
            if ~isempty(obj.Nodesets(1).ID),
                nn=length(obj.Nodesets);
            else
                nn=0;
            end
            fprintf('Number of nodal point sets           = %13d\n',nn)
            if ~isempty(obj.Sidesets(1).ID),
                ns=length(obj.Sidesets);
            else
                ns=0;
            end
            fprintf('Number of element side sets          = %13d\n',ns)
            %%%
            if isempty(obj.GlobalVars(1).Name)
                ng=0;
            else
                ng=length(obj.GlobalVars);
            end
            fprintf('Number of global variables           = %13d\n',ng)
            if isempty(obj.NodalVars(1).Name)
                nn=0;
            else
                nn=length(obj.NodalVars);
            end
            fprintf('Number of variables at each node     = %13d\n',nn)
            if isempty(obj.ElemVars(1,1).Name)
                ne=0;
            else
                ne=length(obj.ElemVars(1,:));
            end
            fprintf('Number of variables at each element  = %13d\n',ne)
            if isempty(obj.NodesetVars(1,1).Name)
                nns=0;
            else
                nns=length(obj.NodesetVars(1,:));
            end
            fprintf('Number of variables at each nodeset  = %13d\n',nns)
            if isempty(obj.SidesetVars(1,1).Name)
                nss=0;
            else
                nss=length(obj.SidesetVars(1,:));
            end
            fprintf('Number of variables at each sideset  = %13d\n\n',nss)
            fprintf('Number of time steps in the database = %8d\n',length(obj.Time))
        end
        function f=limits(obj,blkidx)
            if nargin<2,
                blkidx=1:length(obj.Blocks);
            end
            dim=size(obj.Nodes.Coordinates,2);
            xma=-1e16;
            yma=-1e16;
            zma=-1e16;
            xmi=1e16;
            ymi=1e16;
            zmi=1e16;
            for i=1:length(blkidx),
                xma=max(xma,max(obj.Nodes.Coordinates(obj.Blocks(blkidx(i)).Connectivity,1)));
                xmi=min(xmi,min(obj.Nodes.Coordinates(obj.Blocks(blkidx(i)).Connectivity,1)));
                if dim>1,
                    yma=max(yma,max(obj.Nodes.Coordinates(obj.Blocks(blkidx(i)).Connectivity,2)));
                    ymi=min(ymi,min(obj.Nodes.Coordinates(obj.Blocks(blkidx(i)).Connectivity,2)));
                end
                if dim>2,
                    zma=max(zma,max(obj.Nodes.Coordinates(obj.Blocks(blkidx(i)).Connectivity,3)));
                    zmi=min(zmi,min(obj.Nodes.Coordinates(obj.Blocks(blkidx(i)).Connectivity,3)));
                end
            end
            if nargout, 
                if dim==1,
                    f=[xmi xma xma-xmi];
                elseif dim==2,
                    f=[xmi xma xma-xmi
                        ymi yma yma-ymi];
                else
                    f=[xmi xma xma-xmi
                        ymi yma yma-ymi
                        zmi zma zma-zmi];
                end
            else
                disp('Mesh Limits: ')
                fprintf('Min X = % 12.9e, Max X = % 12.9e, Range = % 12.9e\n',xmi,xma,xma-xmi)
                if dim>1,
                    fprintf('Min Y = % 12.9e, Max Y = % 12.9e, Range = % 12.9e\n',ymi,yma,yma-ymi)
                end
                if dim>2,
                    fprintf('Min Z = % 12.9e, Max Z = % 12.9e, Range = % 12.9e\n',zmi,zma,zma-zmi)
                end
            end
        end
        function v=Volume(obj,blkidx)
            if nargin<2,
                blkidx=1:length(obj.Blocks);
            end
            v=0;
            for i=1:length(blkidx),
                ShpObj=ShapeFactory.CreateShape(obj.Blocks(blkidx(i)).ElementType, ...
                    obj.Blocks(blkidx(i)).Connectivity, ...
                    obj.Nodes.Coordinates);
                v=v+sum(ShpObj.Volume);
            end
        end
        function v=Volume_elem(obj,blkidx)
            if nargin<2,
                blkidx=1:length(obj.Blocks);
            end
            nelem=0;
            for i=1:length(blkidx),
                nelem=nelem+size(obj.Blocks(blkidx(i)).Connectivity,1);
            end
            v=zeros(nelem,1);
            ii=1;
            for i=1:length(blkidx),
                nelem=size(obj.Blocks(blkidx(i)).Connectivity,1);
                ShpObj=ShapeFactory.CreateShape(obj.Blocks(blkidx(i)).ElementType, ...
                    obj.Blocks(blkidx(i)).Connectivity, ...
                    obj.Nodes.Coordinates);
                v(ii:ii+nelem-1)=ShpObj.Volume;
                ii=ii+nelem;
            end
        end
        %%%
        function obj=AddGlobalVar(obj,name,data,time)
            if nargin<4 || isempty(time)
                if ~isempty(obj.Time)
                    time=obj.Time;
                else
                    error('Time has not been yet defined: Exodus.m')
                end
            end
            [timestep]=obj.chkTime(time);
            %%%%%%%%
            if ~isempty(obj.GlobalVars(1).Name)
                for i=1:length(obj.GlobalVars),
                    if strcmp(name,obj.GlobalVars(i).Name),
                        obj.GlobalVars.Data(:,timestep)=data;
                        return
                    end
                end
                nobj=FEMesh.GlobalVars({name});
                nobj.Data=zeros(size(obj.Nodes.Coordinates,1),timestep(end));
                nobj.Data(:,timestep)=data;
                obj.GlobalVars=[obj.GlobalVars nobj];
            else
                obj.GlobalVars(1).Name=name;
                obj.GlobalVars(1).Data=zeros(size(obj.Nodes.Coordinates,1),timestep(end));
                obj.GlobalVars(1).Data(:,timestep)=data;
            end
        end
        function obj=AddNodalVar(obj,name,data,time)
            if nargin<4 || isempty(time)
                if ~isempty(obj.Time)
                    time=obj.Time;
                else
                    error('Time has not been yet defined: Exodus.m')
                end
            end
            [timestep]=obj.chkTime(time);
            %%%%%%%%
            if ~isempty(obj.NodalVars(1).Name)
                for i=1:length(obj.NodalVars),
                    if strcmp(name,obj.NodalVars(i).Name),
                        obj.NodalVars(i).Data(:,timestep)=data;
                        return;
                    end
                end
                nobj=FEMesh.NodalVars({name});
                nobj.Data=zeros(size(obj.Nodes.Coordinates,1),timestep(end));
                nobj.Data(:,timestep)=data;
                obj.NodalVars=[obj.NodalVars nobj];
            else
                obj.NodalVars(1).Name=name;
                obj.NodalVars(1).Data=zeros(size(obj.Nodes.Coordinates,1),timestep(end));
                obj.NodalVars(1).Data(:,timestep(end))=data;
            end
        end
        function obj=AddElemVar(obj,varname,blkidx,data,time)
            if nargin<4 || isempty(time)
                if ~isempty(obj.Time)
                    time=obj.Time;
                else
                    error('Time has not been yet defined: Exodus.m')
                end
            end
            [timestep]=obj.chkTime(time);
            %%%
            nblks=length(obj.Blocks);
            nvars=size(obj.ElemVars,2);
            
            if isempty(obj.ElemVars(1,1).Name),
                names=cell(1,1);
            else
                names=cell(nvars+1,1);
            end
            flag=0;
            for j=1:nvars,
                names{j}=obj.ElemVars(1,j).Name;
                if strcmp(names{j},varname),
                    flag=1;
                    varidx=j;
                    break
                end
            end
            if ~flag
                ids=zeros(nblks,1);
                for i=1:nblks,
                    ids(i)=obj.Blocks(i).ID;
                end
                varidx=length(names);
                nobj=FEMesh.ElemVars(ids,{varname});
                if ~isempty(obj.ElemVars(1,1).Name)
                    obj.ElemVars=[obj.ElemVars nobj];
                else
                    obj.ElemVars=nobj;
                end
            end
            %
            if ~isempty(data),
                obj.ElemVars(blkidx,varidx).Data(:,timestep)=data;
            else
                warning('The Element Variable Data for variable %s in block index %d is empty',varname,blkidx)
            end
        end
        function obj=AddSidesetVar(obj,varname,ssidx,data,time)
            if nargin<4 || isempty(time)
                if ~isempty(obj.Time)
                    time=obj.Time;
                else
                    error('Time has not been yet defined: Exodus.m')
                end
            end
            [timestep]=obj.chkTime(time);
            %%%
            nss=length(obj.Sidesets);
            svars=size(obj.SidesetVars,2);
            
            if isempty(obj.SidesetVars(1,1).Name),
                names=cell(1,1);
            else
                names=cell(nvars+1,1);
            end
            flag=0;
            for j=1:svars,
                names{j}=obj.SidesetVars(1,j).Name;
                if strcmp(names{j},varname),
                    flag=1;
                    varidx=j;
                    break
                end
            end
            if ~flag
                ids=zeros(nss,1);
                for i=1:nss,
                    ids(i)=obj.Sidesets(i).ID;
                end
                varidx=length(names);
                nobj=FEMesh.SidesetVars(ids,{varname});
                if ~isempty(obj.SidesetVars(1,1).Name)
                    obj.SidesetsVars=[obj.SidesetsVars nobj];
                else
                    obj.SidesetVars=nobj;
                end
            end
            %
            if ~isempty(data),
                obj.SidesetVars(ssidx,varidx).Data(:,timestep)=data;
            else
                warning('The Sideset Variable Data for variable %s in sideset index %d is empty',varname,ssidx)
            end
        end
        function evar=getElemVar(fexo,ename,tidx)
            evar=[];
            for i=1:length(fexo.Blocks),
                names=fexo.ElemVars.getNames(i);
                idx=strmatch(ename,names,'exact');
                if isempty(idx),
                    error('variable %s not found',ename);
                end
                evar=cat(1,evar,fexo.ElemVars(i,idx).Data(:,tidx));
            end
        end
        %%%%%
        function nodexyz=ElemCoord(obj,blkidx)
            if nargin<2,
                blkidx=1;
            end
            %% gives the element centroids for a block
            ShpObj=ShapeFactory.CreateShape(obj.Blocks(blkidx).ElementType,obj.Blocks(blkidx).Connectivity(1,:),[]);
            intpt=ShpObj.LocalCentroid;
            x=obj.Nodes.Coordinates(:,1);
            y=obj.Nodes.Coordinates(:,2);
            switch size(obj.Nodes.Coordinates,2),
                case 2
                    z=zeros(size(x));
                case 3
                    z=obj.Nodes.Coordinates(:,3);
            end
            %%
            N=ShpObj.Interpolate(intpt);
            switch size(obj.Nodes.Coordinates,2),
                case 2
                    nodexyz=[x(obj.Blocks(blkidx).Connectivity)*N', ...
                        y(obj.Blocks(blkidx).Connectivity)*N', ...
                        zeros(size(obj.Blocks(blkidx).Connectivity,1),1)];
                    %%%%%%
                case 3
                    nodexyz=[x(obj.Blocks(blkidx).Connectivity)*N', ...
                        y(obj.Blocks(blkidx).Connectivity)*N', ...
                        z(obj.Blocks(blkidx).Connectivity)*N'];
                    %%%%%%
                otherwise
                    error('Unknown Dimension')
            end
        end
        function timestep=chkTime(obj,time)
            if isempty(obj.Time),
                obj.Time=time;
                timestep=1:length(obj.Time);
            else
                timestep=zeros(length(time),1);
                for i=1:length(time)
                    idxhold=find(abs(obj.Time-time(i))<eps*100);
                    if ~isempty(idxhold),
                        if length(idxhold)==1,
                            timestep(i)=idxhold;
                        else
                            error('Multiple similar timesteps are defined')
                        end
                    else
                        obj.Time=[obj.Time(:);time(i)];
                        timestep(i)=length(obj.Time);
                    end
                end
            end
        end
        function obj=RandomizeConn(obj,blkidx)
            if nargin==1,
                blkidx=1:length(obj.Blocks);
            end
            for i=1:length(blkidx),
                Sobj=ShapeFactory.CreateShape(obj.Blocks(blkidx(i)).ElementType, ...
                    obj.Blocks(blkidx(i)).Connectivity,obj.Nodes.Coordinates);
                obj.Blocks(blkidx(i)).Connectivity=Sobj.RandomizeConn(obj.Blocks(blkidx(i)).Connectivity);
            end
        end
        function obj=ReverseElements(obj,blkidx)
            if nargin<2,
                blkidx=1:length(obj.Blocks);
            end
            for i=1:length(blkidx),
                Sobj=ShapeFactory.CreateShape(obj.Blocks(blkidx(i)).ElementType, ...
                    obj.Blocks(blkidx(i)).Connectivity,obj.Nodes.Coordinates);
                IDX=find(Sobj.Volume<0);
                obj.Blocks(blkidx(i)).Connectivity(IDX,:)=Sobj.ReverseElement(obj.Blocks(blkidx(i)).Connectivity(IDX,:));
            end
        end
        function [obj,unodes]=RemoveUnusedNodes(obj)
            %% should check nodesets as well and update them.
            usednodes=[];
            for i=1:length(obj.Blocks),
                if size(obj.Blocks(i).Connectivity,1)==1,
                    usednodes=cat(1,usednodes,unique(obj.Blocks(i).Connectivity)');
                else
                    usednodes=cat(1,usednodes,unique(obj.Blocks(i).Connectivity));
                end
            end
            [usedn]=unique(usednodes);
            a=zeros(length(usedn),1);
            a(usedn)=1:length(usedn);
            unodes=obj.Nodes.Coordinates(~a,:);
            obj.Nodes.Coordinates=obj.Nodes.Coordinates(usedn,:);
            if ~isempty(obj.Nodes.NodeNumMap),
                obj.Nodes.NodeNumMap=obj.Nodes.NodeNumMap(usedn);
            end
            for i=1:length(obj.Blocks),
                if size(obj.Blocks(i).Connectivity,1)==1,
                    obj.Blocks(i).Connectivity=a(obj.Blocks(i).Connectivity)';
                else
                    obj.Blocks(i).Connectivity=a(obj.Blocks(i).Connectivity);
                end
            end
        end
        function [coord]=deformedGeometry(obj,disp_name)
            % coord=obj.deformedGeometry(disp_name);
            %
            %% disp_name would be like Disp so DispX, DispY, and DispZ
            %% would be the displacements of the nodes
            coord=obj.Nodes.Coordinates;
            coord(:,1)=coord(:,1)+obj.NodalVars.getnvar(obj.NodalVars,sprintf('%sX',disp_name));
            if size(coord,2)>1,
                coord(:,2)=coord(:,2)+obj.NodalVars.getnvar(obj.NodalVars,sprintf('%sY',disp_name));
            end
            if size(coord,2)>2,
                coord(:,3)=coord(:,3)+obj.NodalVars.getnvar(obj.NodalVars,sprintf('%sZ',disp_name));
            end
        end
        function [fexo]=Copy(obj)
            %%%
            fexo=FEMesh.Exodus(obj.FloatingPointWordSize,obj.FileSize);
            %%
            fexo.Filename=obj.Filename;
            fexo.Title=obj.Title;
            if ~isempty(obj.QARecords),
                fexo.QARecords=FEMesh.QARecords(obj.QARecords.outputQA);
            else
                fexo.QARecords=FEMesh.QARecords;
            end
            
            fexo.InfoRecords=obj.InfoRecords;
            %%%%
            fexo.Nodes=FEMesh.Nodes(obj.Nodes.Coordinates,obj.Nodes.NodeNumMap,obj.Nodes.Names);
            %
            fexo.CoordinateFrames=obj.CoordinateFrames;
            
            fexo.Blocks=FEMesh.Blocks(obj.Blocks.getBlockIDs);
            for i=1:length(fexo.Blocks),
                fexo.Blocks(i).Connectivity=obj.Blocks(i).Connectivity;
                fexo.Blocks(i).ElementType=obj.Blocks(i).ElementType;
                fexo.Blocks(i).Name=obj.Blocks(i).Name;
                fexo.Blocks(i).Status=obj.Blocks(i).Status;
                fexo.Blocks(i).Attributes=obj.Blocks(i).Attributes;
                fexo.Blocks(i).AttributesName=obj.Blocks(i).AttributesName;
            end
            
            fexo.ElementMap=obj.ElementMap;

            fexo.Nemesis=obj.Nemesis;
            
            fexo.Nodesets=FEMesh.Nodesets(obj.Nodesets.getNodesetIDs);
            if ~isempty(fexo.Nodesets(1).ID),
                for i=1:length(fexo.Nodesets),
                    fexo.Nodesets(i).Nodes=obj.Nodesets(i).Nodes;
                    fexo.Nodesets(i).DistFactors=obj.Nodesets(i).DistFactors;
                    fexo.Nodesets(i).Name=obj.Nodesets(i).Name;
                    fexo.Nodesets(i).Status=obj.Nodesets(i).Status;
                end
            end
            
            fexo.Sidesets=FEMesh.Sidesets(obj.Sidesets.getSidesetIDs);
            if ~isempty(fexo.Sidesets(1).ID),
                for i=1:length(fexo.Sidesets),
                    fexo.Sidesets(i).Elements=obj.Sidesets(i).Elements;
                    fexo.Sidesets(i).Sides=obj.Sidesets(i).Sides;
                    fexo.Sidesets(i).DistFactors=obj.Sidesets(i).DistFactors;
                    fexo.Sidesets(i).Status=obj.Sidesets(i).Status;
                    fexo.Sidesets(i).Name=obj.Sidesets(i).Name;
                end
            end
            
            fexo.Time=obj.Time;
            
            if ~isempty(obj.GlobalVars(1).Name),
                gnames=cell(length(obj.GlobalVars),1);
                for i=1:length(obj.GlobalVars),
                    gnames{i}=obj.GlobalVars(i).Name;
                end
                fexo.GlobalVars=FEMesh.GlobalVars(gnames);
                for i=1:length(gnames),
                    fexo.GlobalVars(i).Data=obj.GlobalVars(i).Data;
                end
            end
            
            if ~isempty(obj.NodalVars(1).Name),
                nnames=cell(length(obj.NodalVars),1);
                for i=1:length(obj.NodalVars),
                    nnames{i}=obj.NodalVars(i).Name;
                end
                fexo.NodalVars=FEMesh.NodalVars(nnames);
                for i=1:length(nnames),
                    fexo.NodalVars(i).Data=obj.NodalVars(i).Data;
                end
            end
            
            if ~isempty(obj.ElemVars(1,1).Name),
                enames=cell(size(obj.ElemVars,2),1);
                ids=zeros(size(obj.ElemVars,1),1);
                for j=1:length(obj.ElemVars(1,:)),
                    enames{j}=obj.ElemVars(1,j).Name;
                end
                for i=1:length(obj.ElemVars(:,1)),
                    ids(i)=obj.ElemVars(i,1).BlockID;
                end
                fexo.ElemVars=FEMesh.ElemVars(ids,enames);
                
                for i=1:length(ids),
                    for j=1:length(enames),
                        fexo.ElemVars(i,j).Data=obj.ElemVars(i,j).Data;
                    end
                end
            end
            
        end
        function [fexo,nelems_per_blk]=MergeBlks(fexo,newblkid,blkidx)
            % this code merges one set of blocks at a time
            if nargin==2,
                blkidx=1:length(fexo.Blocks); % merge all of the blocks
            end
            %%%%
            nelems=0;
            elmtype=fexo.Blocks(blkidx(1)).ElementType;
            nnodes=size(fexo.Blocks(blkidx(1)).Connectivity,2);
            nelems_per_blk=zeros(length(blkidx),1);
            for i=1:length(blkidx),
                [n,m]=size(fexo.Blocks(blkidx(i)).Connectivity);
                if strcmpi(fexo.Blocks(blkidx(i)).ElementType(1:3),elmtype(1:3)) || m~=nnodes,
                    nelems=nelems+n;
                    nelems_per_blk(i)=n;
                else
                    error('Element Types must be the same for merge to work')
                end
            end
            blks=FEMesh.Blocks(newblkid);
            blks.ElementType=elmtype;
            blks.Name=fexo.Blocks(blkidx(1)).Name;
            %% Element Attributes are not mapped as of now
            blks.Connectivity=zeros(nelems,nnodes);
            ii=1;
            for i=1:length(blkidx),
                n=size(fexo.Blocks(blkidx(i)).Connectivity,1);
                blks.Connectivity(ii:ii+n-1,:)=fexo.Blocks(blkidx(i)).Connectivity;
                ii=ii+n;
            end
            
            %
            %  Sidesets are not addressed.
            %
            
            %%% element variables
            if ~isempty(fexo.ElemVars(1,1).Name),
                nms=fexo.ElemVars.getNames(blkidx(1));
                evars=FEMesh.ElemVars(newblkid,nms);
               
                for i=1:length(nms),
                    evars(1,i).Data=zeros(nelems,1);
                end
                for j=1:length(nms),
                    ii=1;
                    for i=1:length(blkidx),
                        n=size(fexo.Blocks(blkidx(i)).Connectivity,1);
                        evars(1,j).Data(ii:ii+n-1)=fexo.ElemVars(blkidx(i),j).Data;
                        ii=ii+n;
                    end
                end
                fexo.ElemVars(blkidx,:)=[];
                fexo.ElemVars(end+1,:)=evars;
            end
            fexo.Blocks(blkidx)=[];
            fexo.Blocks(end+1)=blks;
        end
        function [fexo]=deMergeBlks(fexo,mergedblkidx,unmergedblkid,nelems_per_blk)
            % this code merges one set of blocks at a time
            %%%%
            if sum(nelems_per_blk)~=size(fexo.Blocks(mergedblkidx).Connectivity,1),
                error('nelems_per_blk is not consistent with the given block index')
            end
            %if length(unmergedblkid)==1,
            %    return;
            %end
            blks=FEMesh.Blocks(unmergedblkid);
            ii=1;
            for i=1:length(unmergedblkid),
                blks(i).ElementType=fexo.Blocks(mergedblkidx).ElementType;
                blks(i).Name=fexo.Blocks(mergedblkidx).Name;
                %% Element Attributes are not mapped as of now
                blks(i).Connectivity=fexo.Blocks(mergedblkidx).Connectivity(ii:ii+nelems_per_blk(i)-1,:);
                ii=ii+nelems_per_blk(i);
            end
            
            %
            %  Sidesets are not addressed.
            %
            
            %%% element variables
            if ~isempty(fexo.ElemVars(1,1).Name),
                nms=fexo.ElemVars.getNames;
                evars=FEMesh.ElemVars(unmergedblkid,nms);
                id=fexo.ElemVars.id2idx(fexo.Blocks(mergedblkidx).ID);
                for j=1:length(nms),
                    ii=1;
                    for i=1:length(nelems_per_blk),
                        if ~isempty(fexo.ElemVars(id,j).Data),
                            evars(i,j).Data=fexo.ElemVars(id,j).Data(ii:ii+nelems_per_blk(i)-1);
                            ii=ii+nelems_per_blk(i);
                        end
                    end
                end
                fexo.ElemVars(id,:)=[];
                fexo.ElemVars=[fexo.ElemVars;evars];
            end
            fexo.Blocks(mergedblkidx)=[];
            
            fexo.Blocks=[fexo.Blocks blks];
        end
        function struct = exo_rotate(struct,thetax,thetay,thetaz)
            %EXO_ROTATE
            % Function to rotate Exodus II models in Matlab once loaded using exo_get
            % (T. Simmermacher).  Given rotation angles about the x, y, and z axis, the
            % model is rotated and returned.
            %
            %USAGE:
            %   y = exo_rotate(struct,thetax,thetay,thetaz)
            %
            % where
            %       struct = structure contained Exodus II model (loaded from *.exo
            %               file using exo_get)
            %       thetax = counterclockwise rotation about x-axis, with x-axis point
            %               out of and normal to the clock face,  in radians.
            %       thetay = counterclockwise rotation about y-axis, with y-axis point
            %               out of and normal to the clock face,  in radians.
            %       thetaz = counterclockwise rotation about z-axis, with z-axis point
            %               out of and normal to the clock face,  in radians.
            %            y = rotated structure
            %
            % Created by J. W. Rouse 2/22/10
            
            % INITIALIZE:
            Px = zeros(3);
            Py = Px;
            Pz = Px;
            
            % COUNTERCLOCKWISE ROTATION ABOUT X-AXIS
            %   X-AXIS POINTING OUT OF AND NORMAL TO CLOCK
            Px(1,:) = [1 0 0];
            Px(2,:) = [0 cos(thetax) -sin(thetax)];
            Px(3,:) = [0 sin(thetax) cos(thetax)];
            
            % COUNTERCLOCKWISE ROTATION ABOUT Y-AXIS
            %   Y-AXIS POINTING OUT OF AND NORMAL TO CLOCK
            Py(1,:) = [cos(thetay) 0 sin(thetay)];
            Py(2,:) = [0 1 0];
            Py(3,:) = [-sin(thetay) 0 cos(thetay)];
            
            % COUNTERCLOCKWISE ROTATION ABOUT Z-AXIS
            %   Z-AXIS POINTING OUT OF AND NORMAL TO CLOCK
            Pz(1,:) = [cos(thetaz) -sin(thetaz) 0];
            Pz(2,:) = [sin(thetaz) cos(thetaz) 0];
            Pz(3,:) = [0 0 1];
            
            % ROTATION MATRIX:
            P=Px*Py*Pz;
            
            % ROTATE NODES:
            struct.Nodes.Coordinates = struct.Nodes.Coordinates*P';
            
        end
        function [fexo,ve]=threeD2oneD(obj)
            fexo=FEMesh.Exodus(obj.FloatingPointWordSize,obj.FileSize);
            fexo.Filename=obj.Filename;
            fexo.Title=obj.Title;
            if ~isempty(obj.QARecords),
                fexo.QARecords=FEMesh.QARecords(obj.QARecords.outputQA);
            else
                fexo.QARecords=FEMesh.QARecords;
            end
            
            fexo.InfoRecords=obj.InfoRecords;
            %%%%
            coord=round(obj.Nodes.Coordinates/1e-12)*1e-12;
            
            [nox,iax,ibx]=unique(coord(:,1));
            [noy,iay,iby]=unique(coord(:,2));
            [noz,iaz,ibz]=unique(coord(:,3));
            [ju,dim]=max([length(nox),length(noy),length(noz)]);
            switch dim
                case 1
                    ia=iax;ib=ibx;no=nox;
                    ve='x';
                case 2
                    ia=iay;ib=iby;no=noy;
                    ve='y';
                case 3
                    ia=iaz;ib=ibz;no=noz;
                    ve='z';
            end
            %% set up the nodes
            fexo.Nodes=FEMesh.Nodes(no,[],obj.Nodes.Names);
            %
            fexo.CoordinateFrames=obj.CoordinateFrames;
            %
            % Set up the Blocks
            %
            fexo.Blocks=FEMesh.Blocks(obj.Blocks.getBlockIDs);
            %%
            for i=1:length(obj.ElemVars(:,1)),
                    ids(i)=obj.Blocks(i).ID;
            end
            if ~isempty(obj.ElemVars(1,1).Name),
                enames=cell(size(obj.ElemVars,2),1);
                ids=zeros(size(obj.ElemVars,1),1);
                for j=1:length(obj.ElemVars(1,:)),
                    enames{j}=obj.ElemVars(1,j).Name;
                end
                
                %fexo.ElemVars=FEMesh.ElemVars(ids,enames);
            else
                enames=[];
            end
            
            fexo.ElemVars=FEMesh.ElemVars(ids,enames);
            %%
            for i=1:length(fexo.Blocks),
                %% assume it is a hex
                if ~strcmpi(obj.Blocks(i).ElementType(1:3),'hex'),
                    error('Elements must be all hexs')
                end
                if length(obj.Blocks(i).Connectivity(1,:))~=8,
                    error('Elements must be Hex8s')
                end
                fexo.Blocks(i).ElementType='bar';
                fexo.Blocks(i).Name=obj.Blocks(i).Name;
                fexo.Blocks(i).Status=obj.Blocks(i).Status;
                fexo.Blocks(i).Attributes=[];
                fexo.Blocks(i).AttributesName=[];
                ii=1;
                a=repmat(ia(:),1,8);
                ai=repmat((1:length(ia))',1,8);
                conn=[];
                idx=[];
                for j=1:size(obj.Blocks(i).Connectivity,1),
                    
                    b=repmat(obj.Blocks(i).Connectivity(j,:),length(ia),1);
                    
                    ho=ai(a==b);
                    if ~isempty(ho),
                        conn(ii,:)=ho(:)';
                        idx(ii)=j;
                        ii=ii+1;
                    end
                end
                fexo.Blocks(i).Connectivity=conn;
                for j=1:length(enames),
                    fexo.ElemVars(i,j).Data=obj.ElemVars(i,j).Data(idx,:);
                end
            end
            
            fexo.Nodesets=FEMesh.Nodesets([]);
            fexo.Sidesets=FEMesh.Sidesets([]);
            fexo.Time=obj.Time;
            %
            fexo.Nemesis=[];
            %
            if ~isempty(obj.GlobalVars(1).Name),
                gnames=cell(length(obj.GlobalVars),1);
                for i=1:length(obj.GlobalVars),
                    gnames{i}=obj.GlobalVars(i).Name;
                end
                fexo.GlobalVars=FEMesh.GlobalVars(gnames);
                for i=1:length(gnames),
                    fexo.GlobalVars(i).Data=obj.GlobalVars(i).Data;
                end
            else
                fexo.GlobalVars=FEMesh.GlobalVars([]);
            end
            
            %%
            if ~isempty(obj.NodalVars(1).Name),
                ii=1;
                for i=1:length(obj.NodalVars),
                    nn=obj.NodalVars(i).Name;
                    if strcmpi(nn(end),ve)
                        nnames{ii}=nn(1:end-2);
                        idx(ii)=i;
                        ii=ii+1;
                    end
                end
                fexo.NodalVars=FEMesh.NodalVars(nnames);
                for i=1:length(nnames),
                    fexo.NodalVars(i).Data=obj.NodalVars(idx(i)).Data(ia,:);
                end
            else
                fexo.NodalVars=FEMesh.NodalVars([]);
            end
            %%
            
            
        end
        function [fexo]=oneD2threeD(obj)
            w=.1*ones(size(obj.Nodes.Coordinates,1),1); % used to expand the 1D nodes to 3d
            %
            fexo=FEMesh.Exodus(obj.FloatingPointWordSize,obj.FileSize);
            fexo.Filename=obj.Filename;
            fexo.Title=obj.Title;
            if ~isempty(obj.QARecords),
                fexo.QARecords=FEMesh.QARecords(obj.QARecords.outputQA);
            else
                fexo.QARecords=FEMesh.QARecords;
            end
            
            fexo.InfoRecords=obj.InfoRecords;
            %%
            %% set up the nodes
            
            no=[obj.Nodes.Coordinates -w -w
                obj.Nodes.Coordinates w -w
                obj.Nodes.Coordinates w w
                obj.Nodes.Coordinates -w w];
            nnodes=length(obj.Nodes.Coordinates);
            fexo.Nodes=FEMesh.Nodes(no,[],{'X','Y','Z'});
            fexo.Blocks=FEMesh.Blocks(obj.Blocks.getBlockIDs);
            %%
            ids=zeros(size(obj.ElemVars,1),1);
            for i=1:length(obj.ElemVars(:,1)),
                ids(i)=obj.Blocks(i).ID;
            end
            if ~isempty(obj.ElemVars(1,1).Name),
                enames=cell(size(obj.ElemVars,2),1);
                ids=zeros(size(obj.ElemVars,1),1);
                for j=1:length(obj.ElemVars(1,:)),
                    enames{j}=obj.ElemVars(1,j).Name;
                end
                
            else
                enames=[];
            end
             fexo.ElemVars=FEMesh.ElemVars(ids,enames);
            %%
            for i=1:length(fexo.Blocks),
                %% assume it is a bar element

                fexo.Blocks(i).ElementType='HEX8';
                fexo.Blocks(i).Name=obj.Blocks(i).Name;
                fexo.Blocks(i).Status=obj.Blocks(i).Status;
                fexo.Blocks(i).Attributes=[];
                fexo.Blocks(i).AttributesName=[];
                cn=obj.Blocks(i).Connectivity;
                conn=[cn(:,1) cn(:,1)+nnodes cn(:,1)+2*nnodes cn(:,1)+3*nnodes ...
                    cn(:,2) cn(:,2)+nnodes cn(:,2)+2*nnodes cn(:,2)+3*nnodes];
                fexo.Blocks(i).Connectivity=conn;
                for j=1:length(enames),
                    fexo.ElemVars(i,j).Data=obj.ElemVars(i,j).Data;
                end
            end
            
            fexo.Nodesets=FEMesh.Nodesets([]);
            fexo.Nodesets.ID=1;
            fexo.Nodesets.Nodes=(1:size(fexo.Nodes.Coordinates,1));
            
            fexo.Sidesets=FEMesh.Sidesets([]);
            fexo.Time=obj.Time;
            %
            fexo.Nemesis=[];
            %
            if ~isempty(obj.GlobalVars(1).Name),
                gnames=cell(length(obj.GlobalVars),1);
                for i=1:length(obj.GlobalVars),
                    gnames{i}=obj.GlobalVars(i).Name;
                end
                fexo.GlobalVars=FEMesh.GlobalVars(gnames);
                for i=1:length(gnames),
                    fexo.GlobalVars(i).Data=obj.GlobalVars(i).Data;
                end
            else
                fexo.GlobalVars=FEMesh.GlobalVars([]);
            end
            
            %%
            if ~isempty(obj.NodalVars(1).Name),
                ii=1;
                for i=1:length(obj.NodalVars),
                    nn=obj.NodalVars(i).Name;
                    nnames{ii}=nn(1:end-2);
                    idx(ii)=i;
                    ii=ii+1;
                end
                fexo.NodalVars=FEMesh.NodalVars(nnames);
                for i=1:length(nnames),
                    fexo.NodalVars(i).Data=repmat(obj.NodalVars(idx(i)).Data,4,1);
                end
            else
                fexo.NodalVars=FEMesh.NodalVars([]);
            end
            fexo.NodesetVars=FEMesh.NodesetVars;
            fexo.SidesetVars=FEMesh.SidesetVars;
            %%
            
        end
        function F=plot1d(obj,varname,pas,timestep,dispvar)
            % pas -> pause time (.1)
            %
            % F=plot1d(obj,varname,pas,timestep,dispvar)
            %
            %  F are movie frames
            %
            %   to play =>> movie(figure,F,1,100)
            %
            %
            if nargin<4 || isempty(timestep),
                timestep=1:length(obj.Time);
            end
            if nargin<3 || isempty(pas),
                pas=.01;
            end
            if nargin<5 || isempty(dispvar),
                disp=zeros(size(obj.Nodes.Coordinates,1),length(timestep));
            else
                isnodal=obj.NodalVars.isNodalVar(dispvar);
                if isnodal
                    disp=obj.NodalVars(isnodal).Data(:,timestep);
                end
            end
            node=obj.Nodes.Coordinates(:,1);  % assume it is in the x direction
            isnodal=obj.NodalVars.isNodalVar(varname);
            nblks=length(obj.Blocks);
            cx=zeros(nblks,2);
            name=cell(nblks,1);
            for i=1:nblks,
                cx(i,1)=max(max(obj.Blocks(i).Connectivity));
                cx(i,2)=min(min(obj.Blocks(i).Connectivity));
                name{i}=deblank(obj.Blocks(i).Name');
            end
            cl=zeros(nblks,1,3);
            cl(:,1,:)=obj.color(1:nblks);
            
            int=(.6:(.3/(nblks-1)):.9);
            if isnodal
                nx=-1e16;
                nn=1e16;
                for i=1:length(timestep),
                    nx=max([node+disp(:,i);nx]);
                    nn=min([node+disp(:,i);nn]);
                end
                
                mx=max(max(obj.NodalVars(isnodal).Data));
                mn=min(min(obj.NodalVars(isnodal).Data));
                
                for i=1:length(timestep),
                    nod=node+disp(:,i);
                    
                    if nblks>1,
                        patch(nod([cx cx(:,[2 1])])',repmat([mx mx mn mn]',1,nblks),ones(4,nblks),'FaceColor','flat','CData',cl,'linestyle','none')
                        text(mean(nod(cx),2),int*mx,name,'HorizontalAlignment','center')
                    end
                    alpha(.5)
                    hline=line(nod,obj.NodalVars(isnodal).Data(:,timestep(i)));
                    axis([nn nx mn mx])
                    xlabel('Thickness')
                    ylabel(varname)
                    title(sprintf('Time = %g',obj.Time(timestep(i))))
                    if nargout==1,
                        set(gcf,'Renderer','zbuffer')
                        F(i)=getframe(gcf);
                    end
                    if i<length(timestep),
                        pause(pas)
                        plot(0,0,'w')
                    end
                    
                end
                return
            end
            
            iselem=obj.ElemVars.isElemVar(varname);
            
            if iselem,
                evar=[];
                coord=[];
                cdisp=[];
                for j=1:length(obj.Blocks),
                    conn=obj.Blocks(j).Connectivity;
                    coord=cat(1,coord,(node(conn(:,1))+node(conn(:,2)))/2);
                    cdisp=cat(1,cdisp,(disp(conn(:,1),:)+disp(conn(:,2),:))/2);
                    evar=cat(1,evar,obj.ElemVars(j,iselem).Data(:,timestep));
                end
                nx=-1e16;
                nn=1e16;
                for i=1:length(timestep),
                    nx=max([node+disp(:,i);nx]);
                    nn=min([node+disp(:,i);nn]);
                end
                mx=max(max(evar));
                mn=min(min(evar));
                for i=1:length(timestep),
                    cnod=coord+cdisp(:,i);
                    nod=node+disp(:,i);
                    
                    if nblks>1,
                        patch(nod([cx cx(:,[2 1])])',repmat([mx mx mn mn]',1,nblks),ones(4,nblks),'FaceColor','flat','CData',cl,'linestyle','none')
                        text(mean(nod(cx),2),int*mx,name,'HorizontalAlignment','center')
                    end
                    alpha(.5)
                    hline=line(cnod,evar(:,i));
                    set(hline,'linewidth',2,'color','k');
                    axis([nn nx mn mx])
                    xlabel('Thickness')
                    ylabel(varname)
                    title(sprintf('Time = %g',obj.Time(timestep(i))))
                    if nargout==1,
                        set(gcf,'Renderer','zbuffer')
                        F(i)=getframe(gcf);
                    end
                    if i<length(timestep),
                        pause(pas)
                        plot(0,0,'w')
                    end
                   
                end
                return
            end
            %%%%%%%
            isglobal=obj.GlobalVars.isGlobalVar(varname);
            if isglobal
                plot(obj.Time(timestep),obj.GlobalVars(isglobal).Data)
                xlabel('Time')
                ylabel(varname)
                return
            else
                warning('Variable not found: %s',varname)
            end
        end
        function plotth(obj,varname,num)
            isnodal=obj.NodalVars.isNodalVar(varname);
            if isnodal,
                plot(obj.Time,obj.NodalVars(isnodal).Data(num,:))
                xlabel('Time')
                ylabel(obj.NodalVars(isnodal).Name)
                return
            end
            iselem=obj.ElemVars.isElemVar(varname);
            if iselem,
                [enum,blk]=obj.Blocks.LocalElemNum(num);
                plot(obj.Time,obj.ElemVars(blk,iselem).Data(enum,:))
                xlabel('Time')
                ylabel(obj.ElemVars(blk,iselem).Name)
                return
            end
            isglobal=obj.GlobalVars.isGlobalVar(varname);
            if isglobal
                plot(obj.Time,obj.GlobalVars(isglobal).Data)
                xlabel('Time')
                ylabel(varname)
                return
            else
                warning('Variable not found: %s',varname)
            end
        end
        %%%%
        % The next functions are somewhat flakey
        %%%%
        function h=plot(obj,blknum,color,h)
            if nargin==1 || isempty(blknum),
                blknum=1:length(obj.Blocks);
            end
            if nargin<4,
                h=figure;
                set(h,'Visible','off')
            end
            hold on
            for ii=1:length(blknum)
                i=blknum(ii);
                sobj=ShapeFactory.CreateShape(obj.Blocks(i).ElementType, ...
                    obj.Blocks(i).Connectivity,obj.Nodes.Coordinates);
                if nargin >2,
                    plot(sobj,1:size(obj.Blocks(i).Connectivity,1),color)
                else
                    plot(sobj,1:size(obj.Blocks(i).Connectivity,1),FEMesh.Exodus.color(i));
                end
                
            end
            if size(obj.Nodes.Coordinates,2)==3,
                view(3)
            else
                view(2)
            end
            axis equal
            hold off
            if nargin<4,
                set(h,'Visible','on')
            end
            if nargout==0,clear h;end
        end
        function h=plote(obj,elemnum,color,h)
            [elemnum,blknum]=obj.Blocks.LocalElemNum(elemnum);
            if nargin<4,
                h=figure;
                set(h,'Visible','off')
            end
            hold on
            for ii=1:length(blknum)
                i=blknum(ii);
                sobj=ShapeFactory.CreateShape(obj.Blocks(i).ElementType, ...
                    obj.Blocks(i).Connectivity,obj.Nodes.Coordinates);
                if nargin >2,
                    plot(sobj,elemnum,color)
                else
                    plot(sobj,elemnum,FEMesh.Exodus.color(i));
                end
                
            end
            if size(obj.Nodes.Coordinates,2)==3,
                view(3)
            else
                view(2)
            end
            axis equal
            hold off
            if nargin<4,
                set(h,'Visible','on')
            end
            if nargout==0,clear h;end
        end
        function h=plotns(obj,nsnum,color,h)
            if nargin<4,
                h=figure;
                set(h,'Visible','off')
            end
            if nargin<2 || isempty(nsnum),
                nsnum=1:length(obj.Nodesets);
            end
            hold on
            for ii=1:length(nsnum)
                nsvi=obj.Nodesets(nsnum(ii)).Nodes;
                
                if nargin >2 && ~isempty(color),
                    plot3(obj.Nodes.Coordinates(nsvi,1),obj.Nodes.Coordinates(nsvi,2), ...
                        obj.Nodes.Coordinates(nsvi,3),'*','Color',color(ii))
                else
                    plot3(obj.Nodes.Coordinates(nsvi,1),obj.Nodes.Coordinates(nsvi,2), ...
                        obj.Nodes.Coordinates(nsvi,3),'*','Color',FEMesh.Exodus.color(ii));
                end
                
            end
            if size(obj.Nodes.Coordinates,2)==3,
                view(3)
            else
                view(2)
            end
            axis equal
            hold off
            if nargin<4,
                set(h,'Visible','on')
            end
            if nargout==0,clear h;end
        end
        function fexo=GaussMeshRefinement(obj,blkidxs,numpts,type)
            %
            %
            %  fexo=GaussMeshRefinement(obj,blkidxs,numpts,type)
            %
            %  s=> equal spacing
            %  g=> spacing based upon gauss points
            if nargin<4,
                type='s';
            else
                type=type(1);
            end
            x=obj.Nodes.Coordinates(:,1);
            y=obj.Nodes.Coordinates(:,2);
            if size(obj.Nodes.Coordinates,2)==2,
                z=zeros(size(x));
            else
                z=obj.Nodes.Coordinates(:,3);
            end
            fexo=FEMesh.Exodus;
            fexo.Title='Default Matlab Title';
            fexo.Filename='Matlab.g';
            a=ver('Matlab');
            fexo.QARecords=FEMesh.QARecords('Matlab',sprintf('Version %s%s',a.Version,a.Release));
            fexo.Blocks=FEMesh.Blocks([obj.Blocks(blkidxs).ID]);
            fexo.Nodesets=FEMesh.Nodesets([]);
            fexo.Sidesets=FEMesh.Sidesets([]);
            fexo.GlobalVars=FEMesh.GlobalVars([]);
            fexo.NodalVars=FEMesh.NodalVars([]);
            fexo.ElemVars=FEMesh.ElemVars([obj.Blocks(blkidxs).ID],[]);
            fexo.NodesetVars=FEMesh.NodesetVars([]);
            coordx=[];
            coordy=[];
            coordz=[];
            nelem=0;
            nnode=0;
                        
            for i=1:length(blkidxs),
                Sobj=ShapeFactory.CreateShape(obj.Blocks(blkidxs(i)).ElementType, ...
                    obj.Blocks(blkidxs(i)).Connectivity(1,:),obj.Nodes.Coordinates(1,:));
                [v,conn,elemtype]=Sobj.GenSubCells(numpts(i),type);
                elemtype=obj.Blocks(blkidxs(i)).ElementType;
                N=Sobj.Interpolate(v);
                cx=N*x(obj.Blocks(blkidxs(i)).Connectivity)';
                cy=N*y(obj.Blocks(blkidxs(i)).Connectivity)';
                cz=N*z(obj.Blocks(blkidxs(i)).Connectivity)';
                coordx=cat(1,coordx,cx(:));
                coordy=cat(1,coordy,cy(:));
                coordz=cat(1,coordz,cz(:));
              
                fexo.Blocks(i).ElementType=elemtype;
                fexo.Blocks(i).Connectivity=FEMesh.Exodus.expand_connectivity(conn, ...
                    size(obj.Blocks(blkidxs(i)).Connectivity,1))+nnode;
                nelem=nelem+size(fexo.Blocks(i).Connectivity,1);
                nnode=nnode+numel(cx);
            end
            [coord,idx,jdx]=FEMesh.Exodus.unique_tol([coordx(:) coordy(:) coordz(:)]);
            fexo.Nodes=FEMesh.Nodes(coord,[]);
            for i=1:length(blkidxs),
                fexo.Blocks(i).Connectivity=jdx(fexo.Blocks(i).Connectivity);
            end
            fexo.ElementMap=(1:nelem)';
        end
        function fexo=Equivalence(fexo)
            [fexo.Nodes.Coordinates,idx,jdx]=FEMesh.Exodus.unique_tol(fexo.Nodes.Coordinates);
            if ~isempty(fexo.Nodes.NodeNumMap),
                fexo.Nodes.NodeNumMap=fexo.Nodes.NodeNumMap(idx);
            end
            for i=1:length(fexo.Blocks),
                if size(fexo.Blocks(i).Connectivity,1)==1,
                     fexo.Blocks(i).Connectivity=jdx(fexo.Blocks(i).Connectivity);
                else
                fexo.Blocks(i).Connectivity=jdx(fexo.Blocks(i).Connectivity);
                end
            end
        end
        function fexo=Sideset2Exodus(obj,ssid)
            fexo=FEMesh.Exodus;
            fexo.Title='Default Matlab Title';
            fexo.Filename='Matlab.g';
            fexo.QARecords=FEMesh.QARecords('Matlab',sprintf('Version %s',obj.matversion));
            fexo.Blocks=FEMesh.Blocks(1);
            fexo.Nodes=FEMesh.Nodes(obj.Nodes.Coordinates,[]);
            fexo.Nodesets=FEMesh.Nodesets([]);
            fexo.Sidesets=FEMesh.Sidesets([]);
            fexo.GlobalVars=FEMesh.GlobalVars([]);
            fexo.NodalVars=FEMesh.NodalVars({'DistFactors'});
            fexo.ElemVars=FEMesh.ElemVars([],[]);
            
            idx=obj.Sidesets.id2idx(ssid);
            %
            elem=obj.Sidesets(idx).Elements;
            sides=obj.Sidesets(idx).Sides;
            
            [elemnum,blkidx]=obj.Blocks.LocalElemNum(elem);
            jj=1;
            for i=1:length(obj.Blocks),
                jdx=find(blkidx==i);
                Sobj=ShapeFactory.CreateShape(obj.Blocks(i).ElementType,obj.Blocks(i).Connectivity(1,:),[]);
                lside=zeros(length(jdx),9); %% 9 is the maximum number of nodes in a side possible for exodus
                numnodes=zeros(length(jdx),1);
                for j=1:length(jdx),
                    sid=Sobj.SideDef(sides(jdx(j)));
                    lside(j,1:length(sid))=sid;
                    numnodes(j)=length(sid);
                end
                numnode=unique(numnodes);
                for j=1:length(unique(numnodes)),
                    jidx=find(numnode(j)==numnodes);
                    fexo.Blocks(jj)=FEMesh.Blocks(jj);
                    fexo.Blocks(jj).ElementType=obj.side2element(numnode(j));
                    vi=size(obj.Blocks(jj).Connectivity,1)*(lside(jidx,1:numnode(j))-1)+repmat(elemnum(i==blkidx),1,numnode(j));
                    fexo.Blocks(jj).Connectivity=obj.Blocks(i).Connectivity(vi);
                    jj=jj+1;
                end
            end
            fexo=fexo.RemoveUnusedNodes;
        end
        function fexo=ModelReduce(obj,blkidxs)
            % blkidxs are the blocks to keep
            if length(blkidxs)~=length(unique(blkidxs)),
                error('Some block ids are repeated')
            end
            fexo=FEMesh.Exodus(obj.FloatingPointWordSize,obj.FileSize);
            fexo.Title=sprintf('%s',obj.Title);
            fexo.Filename='junkppp.g';
            a=ver('Matlab');
            qacell=obj.QARecords.outputQA;
            nqa=length(qacell);
            fexo.QARecords=FEMesh.QARecords(qacell);
            fexo.Blocks=FEMesh.Blocks([obj.Blocks(blkidxs).ID]);
            fexo.Time=obj.Time;
            nidx=[];
            bidx=[];
            nelems=0;
            for i=1:length(obj.Blocks),
                nelems=nelems+size(obj.Blocks(i).Connectivity,1);
            end
            for i=1:length(blkidxs),
                fexo.Blocks(i).ElementType=obj.Blocks(blkidxs(i)).ElementType;
                fexo.Blocks(i).Connectivity=obj.Blocks(blkidxs(i)).Connectivity;
                fexo.Blocks(i).Name=obj.Blocks(blkidxs(i)).Name;
                fexo.Blocks(i).Attributes=obj.Blocks(blkidxs(i)).Attributes;
                fexo.Blocks(i).AttributesName=obj.Blocks(blkidxs(i)).AttributesName;
                fexo.Blocks(i).Status=obj.Blocks(blkidxs(i)).Status;
                
                nidx=cat(1,nidx,unique(obj.Blocks(blkidxs(i)).Connectivity(:)));
                bidx=cat(1,bidx,obj.Blocks.GlobalElemNum(blkidxs(i),1:size(obj.Blocks(blkidxs(i)).Connectivity,1))');
            end
            nidx=unique(nidx);
            
            if ~isempty(obj.Nodes.NodeNumMap),
                nnm=obj.Nodes.NodeNumMap(nidx);
            else
                nnm=nidx;
            end
            fexo.Nodes=FEMesh.Nodes(obj.Nodes.Coordinates(nidx,:),nnm,obj.Nodes.Names);
            
            a=zeros(size(obj.Nodes.Coordinates,1),1);
            a(nidx)=1:length(nidx);
            for i=1:length(blkidxs),
                if size(fexo.Blocks(i).Connectivity,1)==1,
                    fexo.Blocks(i).Connectivity=a(fexo.Blocks(i).Connectivity)';
                else
                    fexo.Blocks(i).Connectivity=a(fexo.Blocks(i).Connectivity);
                end
                
            end
            if ~isempty(obj.ElementMap),
                enm=obj.ElementMap(bidx);
            else
                enm=bidx;
            end
            fexo.ElementMap=enm;
            %%%
            if ~isempty(obj.Nodesets(1).ID),
                fexo.Nodesets=FEMesh.Nodesets([obj.Nodesets.ID]);
                for i=1:length(obj.Nodesets),
                    nods=a(obj.Nodesets(i).Nodes);
                    I=find(nods);
                    fexo.Nodesets(i).Nodes=nods(I);
                    if ~isempty(obj.Nodesets(i).DistFactors),
                        fexo.Nodesets(i).DistFactors=obj.Nodesets(i).DistFactors(I);
                    else
                        fexo.Nodesets(i).DistFactors=[];
                    end
                    fexo.Nodesets(i).Status=obj.Nodesets(i).Status;
                    fexo.Nodesets(i).Name=obj.Nodesets(i).Name;
                end
                idx=[];
                for i=1:length(obj.Nodesets)
                    if ~isempty(fexo.Nodesets(i).Nodes)
                        idx=[idx i];
                    end
                end
                if ~isempty(idx),
                    fexo.Nodesets=fexo.Nodesets(idx);
                else
                    fexo.Nodesets=FEMesh.Nodesets;
                end
            else
                fexo.Nodesets=FEMesh.Nodesets;
            end
            %%%%% Sidesets not done
            fexo.Sidesets=FEMesh.Sidesets;
            if ~isempty(obj.Sidesets(1).ID),
                elemlist=zeros(nelems,1);
                elemlist(bidx)=1:length(bidx);
                fexo.Sidesets=FEMesh.Sidesets([obj.Sidesets.ID]);
                
                for i=1:length(obj.Sidesets),
                    fexo.Sidesets(i).Name=obj.Sidesets(i).Name;
                    elemid=obj.Sidesets(i).Elements;
                    idx=find(elemlist(elemid));
                    fexo.Sidesets(i).Elements=elemlist(elemid(idx));
                    fexo.Sidesets(i).Sides=obj.Sidesets(i).Sides(idx);
                    if ~isempty(fexo.Sidesets(i).Elements),
                        [el,bi]=fexo.Blocks.LocalElemNum(fexo.Sidesets(i).Elements);
                        nn=zeros(length(bi),1);
                        
                        for j=1:length(bi),
                            %[node,I1,I2]=side_nodes(elemtype,nodesperelem,side,dim)
                            siden=obj.side_nodes(fexo.Blocks(bi(j)).ElementType,size(fexo.Blocks(bi(j)).Connectivity,2),fexo.Sidesets(i).Sides(j),3);
                            if iscell(siden),
                                if ~isempty(siden{1}),
                                    nn(j)=length(siden{1});
                                else
                                    nn(j)=length(siden{2});
                                end
                            else
                                nn(j)=length(siden);
                            end
                        end
                        fexo.Sidesets(i).DistFactors=ones(sum(nn),1);
                    end
                    
                    fexo.Sidesets(i).Status=obj.Sidesets(i).Status;
                end
                %% remove empty sidesets
                ess=zeros(length(fexo.Sidesets),1);
                for jj=1:length(fexo.Sidesets),
                    if ~isempty(fexo.Sidesets(jj).Elements),
                        ess(jj)=1;
                    end
                end
                
                fexo.Sidesets=fexo.Sidesets(logical(ess));
                if isempty(fexo.Sidesets),
                    fexo.Sidesets=FEMesh.Sidesets;
                end
            else
                fexo.Sidesets=FEMesh.Sidesets;
            end
            %%%%%

            if ~isempty(obj.GlobalVars) && ~isempty(obj.GlobalVars(1).Name),
                names=cell(length(obj.GlobalVars),1);
                for i=1:length(obj.GlobalVars),
                    names{i}=obj.GlobalVars(i).Name;
                end
                fexo.GlobalVars=FEMesh.GlobalVars(names);
                fexo.GlobalVars=FEMesh.GlobalVars(names);
                for i=1:length(obj.GlobalVars),
                    fexo.GlobalVars(i).Data=obj.GlobalVars(i).Data;
                end
            else
                fexo.GlobalVars=FEMesh.GlobalVars;
            end
            
            %%%%%
            if ~isempty(obj.NodalVars) && ~isempty(obj.NodalVars(1).Name),
                names=cell(length(obj.NodalVars),1);
                for i=1:length(obj.NodalVars),
                    names{i}=obj.NodalVars(i).Name;
                end
                fexo.NodalVars=FEMesh.NodalVars(names);
                for i=1:length(obj.NodalVars),
                    fexo.NodalVars(i).Data=obj.NodalVars(i).Data(nidx,:);
                end
            else
                fexo.NodalVars=FEMesh.NodalVars;
                %%%
            end
            %%
            if ~isempty(obj.ElemVars) && ~isempty(obj.ElemVars(1,1).Name),
                names=obj.ElemVars.getNames;
                fexo.ElemVars=FEMesh.ElemVars([obj.Blocks(blkidxs).ID],names);
                for i=1:length(blkidxs),
                    for j=1:length(names),
                        fexo.ElemVars(i,j).Data=obj.ElemVars(blkidxs(i),j).Data;
                    end
                end
                %%%
            else
                fexo.ElemVars=FEMesh.ElemVars;
            end
            fexo.NodesetVars=FEMesh.NodesetVars([]);
            fexo.SidesetVars=FEMesh.SidesetVars([]);
        end
        function obj=Skin(obj,blkidxs,ssid)
            %blkidxs=blkidxs(1);  %% just do one block for now.
            n=zeros(length(blkidxs),1);
            elemtype=obj.Blocks(blkidxs(1)).ElementType;
            nnod=size(obj.Blocks(blkidxs(1)).Connectivity,2);
            for i=1:length(blkidxs),
                if ~strcmpi(elemtype,obj.Blocks(blkidxs(i)).ElementType) || ...
                    size(obj.Blocks(blkidxs(i)).Connectivity,2) ~= nnod,
                    error('All element types must be the same')
                end
                n(i)=size(obj.Blocks(blkidxs(i)).Connectivity,1);
            end
            
            conn=zeros(sum(n),nnod);
            gelemnums=zeros(sum(n),1);
            ii=1;
            for i=1:length(blkidxs),
                conn(ii:(ii+n(i)-1),:)=obj.Blocks(blkidxs(i)).Connectivity;
                gelemnums(ii:(ii+n(i)-1))=obj.Blocks.GlobalElemNum(blkidxs(i),1:obj.Blocks(blkidxs(i)).numelem);
                ii=ii+n(i);
            end

            Sobj=ShapeFactory.CreateShape(obj.Blocks(blkidxs(i)).ElementType, ...
                conn,obj.Nodes.Coordinates);
            [uelems,usides,unperside]=Sobj.UniqueSideNum;
            
            nss=length(obj.Sidesets);
            if ~isempty(obj.Sidesets(1).ID),
                nss=nss+1;
            end
            for i=1:length(obj.Sidesets),
                if ssid==obj.Sidesets(i).ID;
                    error('Sideset is already in use');
                end
            end
            obj.Sidesets(nss)=FEMesh.Sidesets(ssid);
            obj.Sidesets(nss).Elements=gelemnums(uelems);
            obj.Sidesets(nss).Sides=usides;
            obj.Sidesets(nss).DistFactors=ones(sum(unperside),1);
        end
        function obj=Rigid(obj,mass,I)
            nnode=size(obj.Nodes.Coordinates,1);
            if nargin<3 || isempty(I),
                I = eye(3);
            end
            if nargin<2 || isempty(mass),
                mass=1;
            end
            dofmax=6*nnode;
            
            Rm=zeros(dofmax,6);
            %
            Rm(1:6:dofmax,1)=ones(nnode,1)/sqrt(mass);
            Rm(2:6:dofmax,2)=ones(nnode,1)/sqrt(mass);
            Rm(3:6:dofmax,3)=ones(nnode,1)/sqrt(mass);
            %%%
            [V,P]=eig(I);
            noder=obj.Nodes.Coordinates*V;
            [R]=obj.rot_shape(noder);
            %
            R=R/sqrt(P);
            %
            % Rotate Back
            for i=1:round(2*nnode);
                R(3*i-2:3*i,:)=V*R(3*i-2:3*i,:);
            end
            %% fit the rotated modes to the coordinate system modes
            Rm(:,4:6)=R*V';
            obj.Time=zeros(6,1);
            obj.NodalVars=FEMesh.NodalVars({'DispX','DispY','DispZ','RotX','RotY','RotZ'});
            obj.NodalVars(1).Data=Rm(1:6:end,:);
            obj.NodalVars(2).Data=Rm(2:6:end,:);
            obj.NodalVars(3).Data=Rm(3:6:end,:);
            obj.NodalVars(4).Data=Rm(4:6:end,:);
            obj.NodalVars(5).Data=Rm(5:6:end,:);
            obj.GlobalVars=FEMesh.GlobalVars;
            obj.ElemVars=FEMesh.ElemVars;
        end
        function obj=mkHigherOrder(obj,blkids)
            if nargin<2,
                blkids=1:length(obj.Blocks);
            end
            for i=1:length(blkids),
                n=size(obj.Nodes.Coordinates,1);
                sobj=ShapeFactory.CreateShape(obj.Blocks(blkids(i)).ElementType, ...
                    obj.Blocks(blkids(i)).Connectivity,obj.Nodes.Coordinates);
                try
                    [conn,coord]=sobj.mkHigherOrder;
                catch
                    fprintf('Unable to increase order in block %d\n',obj.Blocks(blkids(i)).ID);
                    continue
                end
                obj.Nodes.Coordinates=[obj.Nodes.Coordinates;coord(n+1:end,:)];
                if ~isempty(obj.Nodes.NodeNumMap),
                    m=double(max(obj.Nodes.NodeNumMap));
                    nc=size(coord(n+1:end,:));
                    obj.Nodes.NodeNumMap=[double(obj.Nodes.NodeNumMap);(m+1:m+nc)'];
                end
                obj.Blocks(i).Connectivity=conn;
                obj.Blocks(i).ElementType=obj.Blocks(i).ElementType(1:3);
            end
            obj.NodalVars=FEMesh.NodalVars;
            obj=obj.Equivalence;
        end
        function [area,sidecentroid]=SSArea(fexo,ssidx)
           
            [elemnum,blkidx]=fexo.Blocks.LocalElemNum(fexo.Sidesets(ssidx).Elements);
            ublkidx=unique(blkidx);
            area=zeros(length(fexo.Sidesets(ssidx).Elements),1);
            sidecentroid=zeros(length(fexo.Sidesets(ssidx).Elements),3);
            for i=1:length(ublkidx),
                vi=find(ublkidx(i)==blkidx);
                ShpObj=ShapeFactory.CreateShape(fexo.Blocks(ublkidx(i)).ElementType, ...
                    fexo.Blocks(ublkidx(i)).Connectivity, ...
                    fexo.Nodes.Coordinates);
                [area(vi),sidecentroid(vi,:)]=ShpObj.SurfArea(fexo.Sidesets(ssidx).Sides(vi),elemnum(vi));
            end
        end
        function fexo=SS2NS(fexo,ssidx)
            if nargin<2,
                ssidx=1:length(fexo.Sidesets);
            end
            for i=1:length(ssidx),
                if any([fexo.Nodesets.ID]==fexo.Sidesets(ssidx(i)).ID)
                    warning('Nodeset %d already exists, skipping...',fexo.Sidesets(ssidx(i)).ID);
                end
                [nodes]=fexo.SSNodes(ssidx(i));
                [no,idx]=unique(nodes);
                if isempty(fexo.Nodesets(1).ID),
                    fexo.Nodesets(1).ID=fexo.Sidesets(ssidx(i)).ID;
                    fexo.Nodesets(1).Nodes=no;
                    fexo.Nodesets(1).DistFactors=fexo.Sidesets(ssidx(i)).DistFactors(idx);
                    fexo.Nodesets(1).Name=fexo.Sidesets(ssidx(i)).Name;
                else
                    ns=FEMesh.Nodesets(fexo.Sidesets(ssidx(i)).ID);
                    ns.Nodes=no;
                    ns.DistFactors=fexo.Sidesets(ssidx(i)).DistFactors(idx);
                    ns.Name=fexo.Sidesets(ssidx(i)).Name;
                    fexo.Nodesets=[fexo.Nodesets(:);ns];
                end
            end
        end
        %%%%
        function [blkids,blkidx,nelem]=BlksinSS(fexo,ssidx)
            n=zeros(length(fexo.Blocks),1);
            for i=1:length(fexo.Blocks),
                a=fexo.Blocks.GlobalElemNum(i,1);
                b=fexo.Blocks.GlobalElemNum(i,fexo.Blocks(i).NumElementsBlk);
                n(i)=sum(ismember(fexo.Sidesets(ssidx).Elements,a:b));
            end
            blkidx=find(n);
            blkids=[fexo.Blocks(blkidx).ID];
            nelem=n(blkidx);
        end
        function fexo=mkUniqueSS(fexo,ssidxkeep,ssidxmod)
            for i=1:length(ssidxmod),
                idx1=ssidxkeep;
                idx2=ssidxmod(i);
                [ju,jui,idx]=intersect(fexo.Sidesets(idx1).Elements, ...
                    fexo.Sidesets(idx2).Elements);
                if ~isempty(ju),
                    fexo.Sidesets(idx2).Elements(idx)=[];
                    fexo.Sidesets(idx2).Sides(idx)=[];
                    if ~isempty(fexo.Sidesets(idx2).Elements),
                        [hold,ju,nelem]=fexo.BlksinSS(idx2);
                    else
                        hold=[];
                        ju=[];
                        nelem=0;
                    end
                    ndf=0;
                    for p=1:length(hold),
                        if strcmpi('HEX',fexo.Blocks(ju(p)).ElementType),
                            ndf=ndf+4*nelem(p);
                        elseif strcmpi('TETRA10',fexo.Blocks(ju(p)).ElementType),
                            ndf=ndf+6*nelem(p);
                        end
                    end
                    fexo.Sidesets(idx2).DistFactors=ones(ndf,1);
                end
                    
            end
        end
        function fexo=mkUniqueSS2(fexo,ssidxkeep,ssidxmod)
            for i=1:length(ssidxmod),
                idx1=ssidxkeep;
                idx2=ssidxmod(i);
                Nodes1=fexo.Sidesets(idx1).Nodes;
                Nodes2=fexo.Sidesets(idx2).Nodes;
                numnodes2=fexo.Sidesets(idx2).NumNodes;
                cnode2=cumsum(numnodes2);
                [cnodes,jui,idx]=intersect(Nodes1,Nodes2);
                vihold=[];
                nidx=[];
                if ~isempty(cnodes),
                    for j=1:length(cnodes),
                        nidxl=find(cnodes(j)==Node2);
                        nidx=cat(1,nidx,nidxl);
                        for p=1:length(nidxl),
                            [ju,ij]=sort([nidxl(p);cnode2]);
                            ij=find(ij==1,1,'first');
                            vihold=cat(1,vihold,ij);
                        end
                    end
                    fexo.Sidesets(idx2).Elements(vihold)=[];
                    fexo.Sidesets(idx2).Sides(vihold)=[];
                    numnodes2(vihold)=[];
                    fexo.Sidesets(idx2).DistFactors=ones(sum(numnodes2),1);
                    fexo.Sidesets(idx2).Nodes=[];
                    fexo.Sidesets(idx2).NumNodes=[];
                    [nodes,numnodes]=fexo.SSNodes(idx2);
                    fexo.Sidesets(idx2).Nodes=nodes;
                    fexo.Sidesets(idx2).NumNodes=numnodes;
                end
            end
        end
        function [val,numnodes]=SSNodes(fexo,ssidx)
            val=zeros(length(fexo.Sidesets(ssidx).DistFactors),1);
            numnodes=zeros(length(fexo.Sidesets(ssidx).Elements),1);
            [~,blkidx]=fexo.Blocks.LocalElemNum(fexo.Sidesets(ssidx).Elements);
            dim=size(fexo.Nodes.Coordinates,2);
            nblks=unique(blkidx);
            for i=1:length(nblks),
                data(i).I=find(blkidx==nblks(i));
                [data(i).node,data(i).I1,data(i).I2]=fexo.side_nodes(fexo.Blocks(nblks(i)).ElementType, ...
                    size(fexo.Blocks(nblks(i)).Connectivity,2), ...
                    fexo.Sidesets(ssidx).Sides(data(i).I),dim);
                if ~isempty(data(i).I1),
                    if iscell(data(i).node),
                        numnodes(data(i).I(data(i).I1))=size(data(i).node{1},2);
                    else
                        numnodes(data(i).I(data(i).I1))=size(data(i).node,2);
                    end
                end
                if ~isempty(data(i).I2),
                    numnodes(data(i).I(data(i).I2))=size(data(i).node{2},2);
                end
            end
            cnnodes=[1;cumsum(numnodes)+1];
            for i=1:length(nblks),
                elems=fexo.Sidesets(ssidx).Elements;
                
                [ielems,blkidx]=fexo.Blocks.LocalElemNum(elems(data(i).I));
                
                conn=fexo.Blocks(nblks(i)).Connectivity(ielems,:);
                if ~isempty(data(i).I1),
                    idxstart=cnnodes(data(i).I(data(i).I1));
                    if iscell(data(i).node),
                        nnodes=size(data(i).node{1},2);
                        nodes=data(i).node{1}';
                    else
                        nnodes=size(data(i).node,2);
                        nodes=data(i).node';
                    end
                    
                    idx=repmat(idxstart,1,nnodes)+repmat(0:(nnodes-1),length(idxstart),1);
                    idx=idx';
                    nelem1=length(data(i).I1);
                    cc=repmat(1:nelem1,nnodes,1);
                    vi=nelem1*(nodes(:)-1)+cc(:);
                    val(idx)=conn(vi);
                end
                if ~isempty(data(i).I2),
                    idxstart=cnnodes(data(i).I(data(i).I2));
                    nnodes=size(data(i).node{2},2);
                    idx=repmat(idxstart,1,nnodes)+repmat(1:nnodes,length(idxstart),1)';
                    idx=idx';
                    nelem2=length(data(i).I1);
                    nodes=data(i).nodes{2}';
                    cc=repmat(1:nelem2,nnodes,1);
                    vi=nelem2*(nodes(:)-1)+cc(:);
                    val(idx)=conn(vi);
                end
            end
            val=val(logical(val));
        end
    end
    methods (Static)
        function [node,idx,jdx]=unique_tol(node,tol)
            if nargin<2
                tol=1e-8;
            end
            [node,idx,jdx]=unique(tol*round(node/tol),'rows');
        end
        function [tconn]=expand_connectivity(conn,nelem)
            [n,m]=size(conn);
            tconn=zeros(n*nelem,m);
            mnode=max(max(conn));
            for i=1:n,
                tconn(i:n:end,:)=repmat((0:nelem-1)'*mnode,1,m)+conn(i*ones(nelem,1),:);
            end
           
        end
        function [rgb]=color(i)
            colorrgb=[0 0 0%0 %Black
                0 0 255%1 %Blue
                202 225 255%2 %Gray Blue
                191 239 255%3 %Light Blue
                0 238 238%4 %Cyan
                110 139 61%5 %Dark Olive
                0 100 0%6 %Dark Green
                0 238 0%7 %Green
                255 255 0%8 %Yellow
                255 193 37%9 %Golden Orange
                238 154 0%10 %Orange
                255 0 0%11 %Red
                255 0 255%12 %Magenta
                255 131 250%13 %Light Magenta
                255 181 197%14 %Pink
                255 255 255%15 %White
                218 112 214]/255;%16]; %+User-defined
            vi=[11 8 7 4 1 14 10 5 6 2 12 13 16]+1;
            rgb=colorrgb(vi(mod(i-1,12)+1),:);
        end
        function cstr=color2(i)
            colors=['k';'b';'c';'g';'y';'r';'m'];
            vi=1:7;
            cstr=colors(vi(mod(i-1,7)+1));
        end
        function [node,I1,I2]=side_nodes(elemtype,nodesperelem,side,dim)
         % this is from exgssn.c in the exodus library
         % First convert element to blocknumber
         % then identify an element type with each element
         %
         I1=1:length(side);
         I2=[];
         switch lower(elemtype(1:3))
            case 'cir'
               node=ones(length(side,1),1);
            case 'sph'
               node=ones(length(side,1),1);
            case {'tru','bea'}
               % truss and beams
               n=[1 2 3];
               node=n(side,:);
               if nodesperelem==2,
                  node=node(:,1:2);
               end
            case 'tri'
               % triangle */
               if dim==2,
                   n=[1 2 4
                       2 3 5
                       3 1 6];
                   node=n(side,:);
                  if nodesperelem==3,
                     node=node(:,1:2);
                  end
               else
                  % triangle 3d */ modified for 3 or 6 node only, could add middle node if necessary
                  I1=find(side<3);
                  I2=find(side>2);
                  n1=[1 2 3 4 5 6
                      3 2 1 6 5 4];
                  n2=[1 2 4
                      2 3 5
                      3 1 6];
                  node{1}=n1(side(I1),:);
                  node{2}=n2(side(I2)-2,:);
                  if nodesperelem==3,
                         if ~isempty(I1),
                             node{1}=node{1}(:,1:3);
                         end
                         if ~isempty(I2),
                             node{2}=node{2}(:,1:2);
                         end
                  end
               end
            case 'qua'
               % quad */
               n=[1 2 5
                   2 3 6
                   3 4 7
                   4 1 8];
               node=n(side,:);
               if nodesperelem==4,
                  node=node(:,1:2);
               end
            case 'she'
               % shell */
               n1=[1 2 3 4 5 6 7 8
                   1 4 3 2 8 7 6 5];
               n2=[1 2 5
                   2 3 6
                   3 4 7
                   4 1 8];
               I1=find(side<3);
               I2=find(side>2);
               node{1}=n1(side(I1),:);
               node{2}=n2(side(I2)-2,:);
               if nodesperelem==4,
                  if ~isempty(I1), % element face
                      node{1}=node{1}(:,1:4);
                  end
                  if ~isempty(I2),
                      node{2}=node{2}(:,1:2);
                  end
               end
            case 'tet'
               % tetra */
               n=[1 2 4 5 9 8
                   2 3 4 6 10 9
                   1 4 3 8 10 7
                   1 3 2 7 6 5];
               node=n(side,:);
               if nodesperelem==4,
                  node=node(:,1:3);
               elseif nodesperelem==8,
                  node=node(:,1:4);
               end
            case 'wed'
               % wedge */
               n1=[1 2 5 4 7 11 13 10
                   2 3 6 5 8 12 14 11
                   1 4 6 3 10 15 12 9];
               n2=[1 3 2 9 8 7
                   4 5 6 13 14 15];
               I1=find(side<4);
               I2=find(side>3);
               node{1}=n1(side(I1),:);
               node{2}=n2(side(I2)-3,:);
               if nodesperelem==6,
                  if ~isempty(I1),
                     node=node{1}(:,1:4);
                  end
                  if ~isempty(I2),
                     node=node{2}(:,1:3);
                  end
               end
            case 'hex'
               % hex */
               n=[1 2 6 5 9 14 17 13 26
                   2 3 7 6 10 15 18 14 25
                   3 4 8 7 11 16 19 15 27
                   1 5 8 4 13 20 16 12 24
                   1 4 3 2 12 11 10 9 22
                   5 6 7 8 17 18 19 20 23];
               node=n(side,:);
               if nodesperelem==8,
                  node=node(:,1:4);
               elseif nodesperelem==20,
                  node=node(:,1:8);
               end
            case 'pyr'
               % pyramid */
               n1=[1 2 5 6 11 10
                   2 3 5 7 12 11
                   3 4 5 8 13 12
                   1 5 4 10 13 9];
               n2=[1 4 3 2 9 8 7 6];
               I1=find(side<5);
               I2=find(side==5);
               node{1}=n1(side(I1),:);
               node{2}=n2(side(I2),:);
               if nodesperelem==5
                  if ~isempty(I1),
                     node{1}=node{1}(:,1:3);
                  end
                  if ~isempty(I2),
                     node{2}=node{2}(:,1:4);
                  end
               end
            otherwise
               error('Unknown element, %s',elemtype);
         end
      end
        function type=side2element(nnodes)
            switch nnodes
                case 2
                    type='Bar2';
                case 3
                    type='Tri3';
                case 4
                    type='Quad4';
                case 5
                    type='Unknown';
                case 6
                    type='Tri6';
                case 7
                    type='Unknown';
                case 8
                    type='Quad8';
                case 9
                    type='Unknown';
                otherwise
                    error('Unknown number of sides')
            end
        end
        function [Rm]=rot_shape(nodexyz)
            nnode=size(nodexyz,1);
            dofmax=6*nnode;
            Rm=zeros(dofmax,3);
            Rm(1:6:dofmax,[2 3])=[nodexyz(:,3) -nodexyz(:,2)];
            Rm(2:6:dofmax,[1 3])=[-nodexyz(:,3) nodexyz(:,1)];
            Rm(3:6:dofmax,[1 2])=[nodexyz(:,2) -nodexyz(:,1)];
            
            Rm(4:6:dofmax,1)=ones(nnode,1);
            Rm(5:6:dofmax,2)=ones(nnode,1);
            Rm(6:6:dofmax,3)=ones(nnode,1);
        end
        function [node_elem,nodes]=inv_connect(conn)
            %
            %
            num_nodes_per_element=size(conn,2);
            %
            conn=conn';
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
    end
end
