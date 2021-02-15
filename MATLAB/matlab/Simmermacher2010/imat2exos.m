%
%  Assumes geometry and shape are in global coordinates
%
%  fexo=imat2exo(fem,shp)
%
%
%  once the fexo file is created, use exo_wr to write
%    a matlab "exodus" file
%
%
function fexo=imat2exos(fem,shp)
%
%
if nargin>1,
    shp=xform(shp,fem,[],1);
end
fem=xform(fem,[],1);
[nmap,I]=sort(fem.node.id);

fexo.Attributes.floating_point_word_size=4;
fexo.Attributes.file_size=1;


fexo.QARecords{1}{1}='IMAT';
fexo.QARecords{1}{2}='2.6';
fexo.QARecords{1}{3}=datestr(now,23);
fexo.QARecords{1}{4}=datestr(now,13);

fexo.Nodes.Coordinates=fem.node.coord(I,:);
fexo.Nodes.NodeNumMap=nmap;
fexo.Nodes.Names={'xcoord','ycoord','zcoord'};
fexo.Nodesets=[];
fexo.Sidesets=[];
fexo.ElemVars=[];
fexo.Title='Translated Universal File';
%
ii=1;

if ~isempty(fem.tl),
    for i=1:length(fem.tl.id),
        fexo.Blocks(i).ID=fem.tl.id(i);
        fexo.Blocks(i).Name=fem.tl.desc{i};
        fexo.Blocks(i).ElementType='TRUSS2';
        fexo.Blocks(i).Attributes=[];
        fexo.Blocks(i).AttributesName=[];
        fexo.Blocks(i).Status=1;
        jj=1;
        for j=2:length(fem.tl.conn{i}),
            conn_glob=[fem.tl.conn{i}(j-1) fem.tl.conn{i}(j)];
            if conn_glob(2) && conn_glob(1),
                idx1=find(conn_glob(1)==nmap);
                idx2=find(conn_glob(2)==nmap);
                conn=[idx1 idx2];
                fexo.Blocks(i).Connectivity(jj,:)=conn;
                jj=jj+1;
                ii=ii+1;
            end
        end
        
    end
else
    i=1;
end
if ~isempty(fem.elem.id),
    pp=i;
    eblocks=unique(fem.elem.prop(:,1));
    for i=1:length(eblocks),
        indprop=find(fem.elem.prop(:,1)==eblocks(i));
        diffblks=unique(fem.elem.type(indprop));
        for idiff=1:length(diffblks),
            inddiff=indprop(fem.elem.type(indprop)==diffblks(idiff));
            switch fem.elem.type(inddiff(1))
                case 11 %{'truss2'}
                    name='truss2';
                case 21 %{'bar2'}
                    name='bar2';
                case 91 %{'tri3'}
                    name='tri3';
                case 94 %{'quad4'}
                    name='shell4';
                case 95 %{'quad8'}
                    name='shell8';
                %case {121,122} %{'rigid'}
                %    name='truss2';
                case 115 %{'hex8'}
                    name='hex8';
                case 111 %{'tetra4'}
                    name='tetra4';
                case 112 %{'wedge6'}
                    name='wedge6';
                case 161 %{'sphere'} lumped mass
                    name='sphere';
                otherwise
                    name=[];
            end
            if ~isempty(name),
                if pp==1,
                    fexo.Blocks.ID=fem.elem.prop(inddiff(1));
                elseif idiff==1,
                    fexo.Blocks(pp).ID=fem.elem.prop(inddiff(1));
                else
                    fexo.Blocks(pp).ID=str2num(sprintf('%d000%d',idiff,fem.elem.prop(inddiff(1))));
                end
                fexo.Blocks(pp).Name=[];
                fexo.Blocks(pp).ElementType=name;
                fexo.Blocks(pp).Attributes=[];
                fexo.Blocks(pp).AttributesName=[];
                fexo.Blocks(pp).Status=1;
                fexo.Blocks(pp).Connectivity=zeros(length(inddiff),length(fem.elem.conn{inddiff(1)}));
                conn=zeros(1,length(fem.elem.conn{inddiff(1)}));
                for j=1:length(inddiff),
                    conn=conn*0;
                    for k=1:length(conn),
                        conn(k)=find(fem.elem.conn{inddiff(j)}(k)==fem.node.id);
                    end
                    fexo.Blocks(pp).Connectivity(j,:)=conn;
                    ii=ii+1;
                end
                pp=pp+1;
            else
                fprintf('Ideas Element Type %d not converted\n',fem.elem.type(inddiff(1)))
            end
        end
    end
end

%% need to update element map to reflect fem.elem.id
fexo.ElementMap=(1:(ii-1))';%% just added
if nargin==2,
    %%%%
    fexo.GlobalVars(1).Name='ModeNumber';
    fexo.GlobalVars(2).Name='EigenFrequency';
    fexo.GlobalVars(2).Name='EigenVectScale';
    
    fexo.GlobalVars(1).Data=(1:length(shp));
    fexo.GlobalVars(2).Data=shp.Frequency';
    
    max_coord=max(max(fexo.Nodes.Coordinates)-min(fexo.Nodes.Coordinates));
    [b,I]=max(abs(shp.Shape));
    b=diag(shp.Shape(I,:))';
    scale=max_coord./b/10;
    fexo.GlobalVars(3).Data=scale;
    for i=1:length(shp),
        shp.Shape(:,i)=shp.Shape(:,i)*scale(i);
    end
    %%%%
    datax=shp.Shape(1:3:end,:);
    datay=shp.Shape(2:3:end,:);
    dataz=shp.Shape(3:3:end,:);
    [common_nodes,I1]=intersect(shp.Node(:,1),nmap);
    [ju,I2]=intersect(nmap,common_nodes);
    
    fexo.Time=shp.Frequency;
    nmodes=length(fexo.Time);
    fexo.NodalVars(1).Name='DispX';
    fexo.NodalVars(2).Name='DispY';
    fexo.NodalVars(3).Name='DispZ';
    
    fexo.NodalVars(1).Data=zeros(length(nmap),nmodes);
    fexo.NodalVars(1).Data(I2,:)=datax(I1,:);
    
    
    fexo.NodalVars(2).Data=zeros(length(nmap),nmodes);
    fexo.NodalVars(2).Data(I2,:)=datay(I1,:);
    
    fexo.NodalVars(3).Data=zeros(length(nmap),nmodes);
    fexo.NodalVars(3).Data(I2,:)=dataz(I1,:);
    
else
    fexo.GlobalVars=[];

    fexo.NodalVars=[];
    fexo.Time=[];
end

