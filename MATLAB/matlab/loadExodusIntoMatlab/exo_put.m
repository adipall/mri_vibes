function exo_put(filename,fexo)
%     
%  exo_put(filename,exo_struct);
%
%  where filename is a string providing the filename of interst
%        exo_struct is a matlab structure from exo_get or equivalent
%


% Author Todd Simmermacher
    % Tel: 844-0096
    % Last Updated: 04/14/2011 - Todd Simmermacher
    % Sandia National Laboratories
    %
   
%

%% Estimate number of dimensions required
ndims=12;  %% upper estimate of nominal number of dims
ndims=ndims+3*length(fexo.Blocks);
if isfield(fexo,'Sidesets') && ~isempty(fexo.Sidesets),
    ndims=ndims+2*length(fexo.Sidesets);
end
if isfield(fexo,'Nodesets') && ~isempty(fexo.Nodesets),
    ndims=ndims+length(fexo.Nodesets);
end
if ndims>=netcdf.getConstant('NC_MAX_DIMS'),
    warning('Model may have too many features for Matlab \nto write the Exodus file: num_dims approximately %d',ndims)
end

obj=ExodusIO(filename,'create');

if isfield(fexo,'Nodesets') && ~isempty(fexo.Nodesets) ...
        || (isa(fexo,'FEMesh.Exodus') && ~isempty(fexo.Nodesets(1).ID)),
    nns=length(fexo.Nodesets);
else
     nns=0;
end

if isfield(fexo,'Sidesets') && ~isempty(fexo.Sidesets) ...
        || (isa(fexo,'FEMesh.Exodus') && ~isempty(fexo.Sidesets(1).ID)),
    sss=length(fexo.Sidesets);
else
    sss=0;
end

num_elem=0;

for i=1:length(fexo.Blocks),
    num_elem=num_elem+size(fexo.Blocks(i).Connectivity,1);
end

if isfield(fexo,'Attributes'),
    obj=obj.create(fexo.Attributes.floating_point_word_size,fexo.Attributes.file_size,...
        fexo.Attributes.last_time_written,fexo.Title,size(fexo.Nodes.Coordinates,2), ...
        size(fexo.Nodes.Coordinates,1),num_elem,length(fexo.Blocks),nns,sss);
else
    obj=obj.create(fexo.FloatingPointWordSize,fexo.FileSize,fexo.LastTimeWritten, ...
        fexo.Title,size(fexo.Nodes.Coordinates,2),size(fexo.Nodes.Coordinates,1), ...
        num_elem,length(fexo.Blocks), ...
        nns,sss);
end

obj=obj.wr_Time(fexo.Time);


if isfield(fexo,'Nemesis'),
    obj=obj.wr_Nemesis(fexo.Nemesis);
end

         
if isfield(fexo,'InfoRecords') && ~isempty(fexo.InfoRecords),
    obj=obj.wr_Info(fexo.InfoRecords);
end

if isfield(fexo,'QARecords') && ~isempty(fexo.QARecords),
    obj=obj.wr_QA(fexo.QARecords);
end

obj=obj.wr_Nodes(fexo.Nodes.Coordinates,fexo.Nodes.NodeNumMap,fexo.Nodes.Names);

% 
nb=length(fexo.Blocks);
ids=zeros(nb,1);
conn=cell(nb,1);
elemtype=cell(nb,1);
blkname=cell(nb,1);
status=zeros(nb,1);
att=cell(nb,1);
attname=cell(nb,1);

for i=1:nb,
    ids(i)=fexo.Blocks(i).ID;
    conn{i}=fexo.Blocks(i).Connectivity;
    elemtype{i}=fexo.Blocks(i).ElementType;
    blkname{i}=fexo.Blocks(i).Name;
    status(i)=fexo.Blocks(i).Status;
    att{i}=fexo.Blocks(i).Attributes;
    attname{i}=fexo.Blocks(i).AttributesName;
end

obj=obj.wr_Blocks(ids,conn,elemtype,blkname,status,att,attname);

if isfield(fexo,'ElementMap') && ~isempty(fexo.ElementMap) ...
    || (isa(fexo,'FEMesh.Exodus') && ~isempty(fexo.ElementMap)),
    obj=obj.wr_EMap(fexo.ElementMap);
end

if isfield(fexo,'Nodesets') && ~isempty(fexo.Nodesets) ...
    || (isa(fexo,'FEMesh.Exodus') && ~isempty(fexo.Nodesets(1).ID)),
    nb=length(fexo.Nodesets);
    ids=zeros(nb,1);
    nodes=cell(nb,1);
    status=zeros(nb,1);
    name=cell(nb,1);
    dist_fact=cell(nb,1);
    
    for i=1:nb,
        ids(i)=fexo.Nodesets(i).ID;
        nodes{i}=fexo.Nodesets(i).Nodes;
        name{i}=fexo.Nodesets(i).Name;
        status(i)=fexo.Nodesets(i).Status;
        dist_fact{i}=fexo.Nodesets(i).DistFactors;
    end
    
    obj=obj.wr_Nodesets(ids,nodes,dist_fact,name,status);
end

if isfield(fexo,'Sidesets') && ~isempty(fexo.Sidesets) ...
    || (isa(fexo,'FEMesh.Exodus') && ~isempty(fexo.Sidesets(1).ID)),
    nb=length(fexo.Sidesets);
    ids=zeros(nb,1);
    elems=cell(nb,1);
    sides=cell(nb,1);
    status=zeros(nb,1);
    name=cell(nb,1);
    dist_fact=cell(nb,1);
    
    for i=1:nb,
        ids(i)=fexo.Sidesets(i).ID;
        elems{i}=fexo.Sidesets(i).Elements;
        sides{i}=fexo.Sidesets(i).Sides;
        name{i}=fexo.Sidesets(i).Name;
        status(i)=fexo.Sidesets(i).Status;
        dist_fact{i}=fexo.Sidesets(i).DistFactors;
    end
    
    obj=obj.wr_Sidesets(ids,elems,sides,dist_fact,name,status);
end

%% write superelement stuff if necessary
if isfield(fexo,'SuperElement') ...
        || (isa(fexo,'FEMesh.Exodus') && ~isempty(fexo.SuperElement)),
    if ~isempty(fexo.SuperElement)  && ~isempty(fexo.SuperElement.CBMap),
        obj=obj.wr_SuperElement(fexo.SuperElement.Kr,fexo.SuperElement.Cr, ...
            fexo.SuperElement.Mr,fexo.SuperElement.CBMap,fexo.SuperElement.OTM, ...
            fexo.SuperElement.OutMap,fexo.SuperElement.OTME,fexo.SuperElement.OutElemMap, ...
            fexo.SuperElement.NumInterfaceNodes,fexo.SuperElement.NumEig, ...
            fexo.SuperElement.NumConstraints);
    end
end

if isfield(fexo,'GlobalVars') && ~isempty(fexo.GlobalVars) ...
        || (isa(fexo,'FEMesh.Exodus') && ~isempty(fexo.GlobalVars(1).Name)),
    if isempty(fexo.Time),
        error('Time vector needs to be defined')
    end
    names=cell(length(fexo.GlobalVars),1);
    data=zeros(length(fexo.Time),length(fexo.GlobalVars));
    for i=1:length(fexo.GlobalVars),
        names{i}=fexo.GlobalVars(i).Name;
        data(:,i)=fexo.GlobalVars(i).Data(:);
    end
    obj=obj.init_GlobalVars(names);
    
    obj=obj.wr_GlobalVars(data');
end

if isfield(fexo,'NodalVars') && ~isempty(fexo.NodalVars) ...
        || (isa(fexo,'FEMesh.Exodus') && ~isempty(fexo.NodalVars(1).Name)),
    if isempty(fexo.Time),
        error('Time vector needs to be defined')
    end
    names=cell(length(fexo.NodalVars),1);
    data=zeros(size(fexo.Nodes.Coordinates,1),length(fexo.NodalVars),length(fexo.Time));
   
    for i=1:length(fexo.NodalVars),
        names{i}=fexo.NodalVars(i).Name;
        data(:,i,:)=fexo.NodalVars(i).Data;
    end
    obj=obj.init_NodalVars(names);
    obj=obj.wr_NodalVars(data);
end

if isa(fexo,'FEMesh.Exodus') && ~isempty(fexo.ElemVars(1,1).Name),
    if isempty(fexo.Time),
        error('Time vector needs to be defined')
    end
    nvars=length(fexo.ElemVars(1,:));
    nblks=length(fexo.Blocks);
    names=cell(nvars,1);
    data=cell(nvars,nblks);
    tt=zeros(nvars,nblks);
    for i=1:nblks,
        for j=1:nvars,
            if ~isempty(fexo.ElemVars(i,j).Data),
                tt(j,i)=1;
            end
            data{j,i}=fexo.ElemVars(i,j).Data;
            if i==1,
                names{j}=fexo.ElemVars(i,j).Name;
            end
        end
    end
    obj=obj.init_ElemVars(names,tt);
    obj=obj.wr_ElemVars(data);
    
elseif (isfield(fexo,'ElemVars') && ~isempty(fexo.ElemVars))
    if isempty(fexo.Time),
        error('Time vector needs to be defined')
    end
    nvars=length(fexo.ElemVars(1).ElemVarData);
    nblks=length(fexo.Blocks);
    names=cell(nvars,1);
    data=cell(nvars,nblks);
    tt=zeros(nvars,nblks);
    for i=1:nblks,
        for j=1:nvars,
            if ~isempty(fexo.ElemVars(i).ElemVarData(j).Data),
                tt(j,i)=1;
            end
            data{j,i}=fexo.ElemVars(i).ElemVarData(j).Data;
            if i==1,
                names{j}=fexo.ElemVars(i).ElemVarData(j).Name;
            end
        end
    end
    obj=obj.init_ElemVars(names,tt);
    obj=obj.wr_ElemVars(data);
    
end

if (isa(fexo,'FEMesh.Exodus') && ~isempty(fexo.NodesetVars(1,1).Name)),
    if isempty(fexo.Time),
        error('Time vector needs to be defined')
    end
    nvars=length(fexo.NodesetVars(1));
    nns=length(fexo.Nodesets);
    names=cell(nvars,1);
    data=cell(nvars,nns);
    tt=zeros(nvars,nns);
    for i=1:nns,
        for j=1:nvars,
            if ~isempty(fexo.NodesetVars(i,j).Data),
                tt(j,i)=1;
            end
            data{j,i}=fexo.NodesetVars(i,j).Data;
            if i==1,
                names{j}=fexo.NodesetVars(i,j).Name;
            end
        end
    end
    obj=obj.init_NodesetVars(names,tt);
    obj=obj.wr_NodesetVars(data);
elseif isfield(fexo,'NodesetVars') && ~isempty(fexo.NodesetVars),
    if isempty(fexo.Time),
        error('Time vector needs to be defined')
    end
    nvars=length(fexo.NodesetVars(1).NodesetVarData);
    nns=length(fexo.Nodesets);
    names=cell(nvars,1);
    data=cell(nvars,nns);
    tt=zeros(nvars,nns);
    for i=1:nns,
        for j=1:nvars,
            if ~isempty(fexo.NodesetVars(i).NodesetVarData(j).Data),
                tt(j,i)=1;
            end
            data{j,i}=fexo.NodesetVars(i).NodesetVarData(j).Data;
            if i==1,
                names{j}=fexo.NodesetVars(i).NodesetVarData(j).Name;
            end
        end
    end
    obj=obj.init_NodesetVars(names,tt);
    obj=obj.wr_NodesetVars(data);
end

% -------------------------------------------------------------------------
% Mike Added
if (isa(fexo,'FEMesh.Exodus') && ~isempty(fexo.SidesetVars(1,1).Name)),
    if isempty(fexo.Time),
        error('Time vector needs to be defined')
    end
    svars=size(fexo.SidesetVars,2);
    nss=length(fexo.Sidesets);
    names=cell(svars,1);
    data=cell(svars,nss);
    tt=zeros(svars,nss);
    for i=1:nss,
        for j=1:svars,
            if ~isempty(fexo.SidesetVars(i,j).Data),
                tt(j,i)=1;
            end
            data{j,i}=fexo.SidesetVars(i,j).Data;
            if i==1,
                names{j}=fexo.SidesetVars(i,j).Name;
            end
        end
    end
    obj=obj.init_SidesetVars(names,tt);
    obj=obj.wr_SidesetVars(data);
    
elseif isfield(fexo,'SidesetVars') && ~isempty(fexo.SidesetVars),
    if isempty(fexo.Time),
        error('Time vector needs to be defined')
    end
    svars=length(fexo.SidesetVars(1).SidesetVarData);
    nss=length(fexo.Sidesets);
    names=cell(svars,1);
    data=cell(svars,nss);
    tt=zeros(svars,nss);
    for i=1:nss,
        for j=1:svars,
            if ~isempty(fexo.SidesetVars(i).SidesetVarData(j).Data),
                tt(j,i)=1;
            end
            data{j,i}=fexo.SidesetVars(i).SidesetVarData(j).Data;
            if i==1,
                names{j}=fexo.SidesetVars(i).SidesetVarData(j).Name;
            end
        end
    end
    obj=obj.init_SidesetVars(names,tt);
    obj=obj.wr_SidesetVars(data);
    
end

% 
obj.close;
