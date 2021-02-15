%
%  for now makes bar2 tracelines
%
%
function fem=exo2imats(filename)

%
color=[0 %Black
       1 %Blue
       2 %Gray Blue 
       3 %Light Blue
       4 %Cyan
       5 %Dark Olive
       6 %Dark Green
       7 %Green 
       8 %Yellow
       9 %Golden Orange
       10 %Orange
       11 %Red
       12 %Magenta
       13 %Light Magenta
       14 %Pink
       15 %White
       16]; %+User-defined
vi=[11 8 7 4 1 14 10 5 6 2 12 13]+1;
color=color(vi); % sort colors so that they correspond approximately to ensight
%
if ischar(filename),
    fexo=exo_get(filename);
else
    fexo=filename;
end

fem.cs.part={8,'Part1'};
fem.cs.id=1;
fem.cs.name={'CS1'};
fem.cs.color=8;
fem.cs.type=0;
fem.cs.matrix=[1 0 0;0 1 0;0 0 1;0 0 0];
%%
if ~isempty(fexo.Nodes.NodeNumMap),
    fem.node.id=fexo.Nodes.NodeNumMap;
else
    fem.node.id=(1:size(fexo.Nodes.Coordinates,1))';
end
fem.node.cs=ones(size(fem.node.id,1),3);
fem.node.color=1*ones(size(fem.node.id,1),1);
fem.node.coord=fexo.Nodes.Coordinates;
%%
ii=1;
jj=1;
for i=1:length(fexo.Blocks),
    for j=1:size(fexo.Blocks(i).Connectivity,1),
        fem.elem.id(jj)=j;
        fem.elem.color(jj)=color(ii);
        if ~isempty(fexo.Nodes.NodeNumMap),
            fem.elem.conn{jj}=fem.node.id(fexo.Blocks(i).Connectivity(j,:));
        else
            fem.elem.conn{jj}=fexo.Blocks(i).Connectivity(j,:);
        end
        switch lower(fexo.Blocks(i).ElementType)
            case {'truss2'}
                fem.elem.type(jj)=11;
            case {'bar2'}
                fem.elem.type(jj)=21;
            case {'tri3'}
                fem.elem.type(jj)=91;
            case {'quad4'}
                fem.elem.type(jj)=94;
            case {'quad8'}
                fem.elem.type(jj)=95;
            case {'rigid'}
                fem.elem.type(jj)=121; %% rigid bar
                fem.elem.type(jj)=122; %% rigid element
            case {'hex8'}
                fem.elem.type(jj)=115;
            case {'tetra4'}
                fem.elem.type(jj)=111;
            case {'wedge6'}
                fem.elem.type(jj)=112;
            case {'sphere'}
                fem.elem.type(jj)=161;  % lumped mass
        end
        fem.elem.prop(jj,:)=[1 1];
        jj=jj+1;
        if ii>length(color),
           ii=1;
        end
    end
    ii=ii+1;
end
fem.elem.type=fem.elem.type';
fem.elem.id=fem.elem.id';
fem.elem.color=fem.elem.color';
fem.elem.conn=fem.elem.conn';
fem.elem.beamdata=zeros(length(fem.elem.id),3);
fem.tl.id=[];
fem.tl.color=[];
fem.tl.desc=[];
fem.tl.conn=[];
%%
[fem.node.id,I]=sort(fem.node.id);
fem.node.coord=fem.node.coord(I,:);

