function [fexo]=exo_get(filename)
    %
    %  fexo=exo_get;
    %
    %  fexo=exo_get(filename);
    %
    % where fexo is a exodus structure and filename is a string
    %
    %
    
    % Author Todd Simmermacher
    % Tel: 844-0096
    % Last Updated: 07/22/2010
    % Sandia National Laboratories
    %
    
    if nargin<1 || isempty(filename);
        [fname, fpath] = uigetfile('*','MultiSelect','off');
        filename=fullfile(fpath,fname);
        
        if isequal(fname,0) || isequal(fpath,0)
            fexo=-1;
            return
        end
    else
        [fpath,fname,fext] = fileparts(filename);
        if isempty(fpath);
            fpath = pwd;
        end
        filename=fullfile(fpath,sprintf('%s%s',fname,fext));
    end
    
    
    obj=ExodusIO(filename);
    
    fexo.Title=obj.rd_Title;
    [fws,fs,ver,api_ver,ltw]=obj.rd_FileAttr;
    fexo.Attributes.floating_point_word_size=fws;
    fexo.Attributes.file_size=fs;
    fexo.Attributes.last_time_written=ltw;  % used for restarts.  If empty, then file isn't a restart file
    fexo.QARecords=obj.rd_QA;
    fexo.InfoRecords=obj.rd_Info;
    [fexo.Nodes.Coordinates,fexo.Nodes.NodeNumMap,fexo.Nodes.Names]=obj.rd_Nodes;
    fexo.CoordFrames=obj.rd_CoordFrames;
    
    ids=obj.rd_BlockIDs;
    for i=1:length(ids),
        fexo.Blocks(i).ID=ids(i);
        [fexo.Blocks(i).Connectivity,fexo.Blocks(i).ElementType,fexo.Blocks(i).Name,fexo.Blocks(i).Status,fexo.Blocks(i).Attributes,fexo.Blocks(i).AttributesName]=obj.rd_Block(ids(i));
        fexo.Blocks(i).NumElementsBlk=size(fexo.Blocks(i).Connectivity,1);
    end
    
    fexo.ElementMap=obj.rd_EMap;
    
    nem=obj.rd_Nemesis;
    if ~isempty(nem),
        fexo.Nemesis=nem;
    end
    
    ids=obj.rd_NodesetIDs;
    for i=1:length(ids),
        fexo.Nodesets(i).ID=ids(i);
        [fexo.Nodesets(i).Nodes,fexo.Nodesets(i).DistFactors,fexo.Nodesets(i).Name,fexo.Nodesets(i).Status]=obj.rd_Nodeset(ids(i));
    end
    
    ids=obj.rd_SidesetIDs;
    for i=1:length(ids),
        fexo.Sidesets(i).ID=ids(i);
        [fexo.Sidesets(i).Elements,fexo.Sidesets(i).Sides,fexo.Sidesets(i).DistFactors,fexo.Sidesets(i).Status,fexo.Sidesets(i).Name]=obj.rd_Sideset(ids(i));
    end
    
    fexo.Time=obj.rd_Time;
    
    %% read superelements if necessary
    [kr,cr,mr,cb,otm,om,otme,oem,ddof,din,deig,dnc]=obj.rd_SuperElement;
    if ~isempty(cb),
        fexo.SuperElement.Kr=kr;
        fexo.SuperElement.Cr=cr;
        fexo.SuperElement.Mr=mr;
        fexo.SuperElement.CBMap=cb;
        fexo.SuperElement.OTM=otm;
        fexo.SuperElement.OutMap=om;
        fexo.SuperElement.OTME=otme;
        fexo.SuperElement.OutElemMap=oem;
        fexo.SuperElement.NumDOF=ddof;
        fexo.SuperElement.NumInterfaceNodes=din;
        fexo.SuperElement.NumEig=deig;
        fexo.SuperElement.NumConstraints=dnc;
    end
    %%
    gnames=obj.rd_GlobalVarsNames;
    if ~isempty(gnames),
        for i=1:length(gnames),
            fexo.GlobalVars(i).Name=gnames{i};
            fexo.GlobalVars(i).Data=obj.rd_GlobalVar(gnames{i});
        end
    end
    
    nnames=obj.rd_NodalVarsNames;
    if ~isempty(nnames),
        for i=1:length(nnames),
            fexo.NodalVars(i).Name=nnames{i};
            fexo.NodalVars(i).Data=obj.rd_NodalVar(nnames{i});
        end
    end
    
    enames=obj.rd_ElemVarsNames;
    if ~isempty(enames),
        for i=1:length(fexo.Blocks),
            fexo.ElemVars(i).ID=fexo.Blocks(i).ID;
            for j=1:length(enames),
                fexo.ElemVars(i).ElemVarData(j).Name=enames{j};
                fexo.ElemVars(i).ElemVarData(j).Data=obj.rd_ElemVar(enames{j},fexo.Blocks(i).ID);
            end
        end
    end
    
    nsnames=obj.rd_NodesetVarsNames;
    if ~isempty(nsnames),
        for i=1:length(fexo.Nodesets),
            fexo.NodesetVars(i).ID=fexo.Nodesets(i).ID;
            for j=1:length(nsnames),
                fexo.NodesetVars(i).NodesetVarData(j).Name=nsnames{j};
                fexo.NodesetVars(i).NodesetVarData(j).Data=obj.rd_NodesetVar(nsnames{j},fexo.Nodesets(i).ID);
            end
        end
    end
    
    %----------------------------------------------------------------------
    % Mike Added
    ssnames=obj.rd_SidesetVarsNames;
    if ~isempty(ssnames),
        for i=1:length(fexo.Sidesets),
            fexo.SidesetVars(i).ID=fexo.Sidesets(i).ID;
            for j=1:length(ssnames),
                fexo.SidesetVars(i).SidesetVarData(j).Name=ssnames{j};
                fexo.SidesetVars(i).SidesetVarData(j).Data=obj.rd_SidesetVar(ssnames{j},fexo.Sidesets(i).ID);
            end
        end
    end
    
    obj.close
end
