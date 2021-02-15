%
% flag = 0 => everything gets read
% flag = 1 => elemvars don't get read
% flag = 2 => elemvars and nodal vars don't get read
% flag = 3 => elemvars, nodalvars, and global vars don't get read
% flag = 4 => elemvars, nodalvars, global vars, and sidesets don't get read
% flag = 5 => elemvars, nodalvars, global vars, sidesets, and nodesets don't get read
% flag = 6 => elemvars, nodalvars, global vars, sidesets, and nodesets don't get read, time info does get read
%
%[fexo]=exo_rd(filename,flag)
function [fexo]=exo_rd(filename,flag)
    if nargin<2,
        flag=0;
    end
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

%%%
[fws,fs,ver,api_ver,ltw]=obj.rd_FileAttr;
fexo=FEMesh.Exodus(fws,fs);
fexo.LastTimeWritten=ltw;
%%
fexo.Filename=fname;
fexo.Title=obj.rd_Title;
if ~isempty(obj.rd_QA),
    fexo.QARecords=FEMesh.QARecords(obj.rd_QA);
else
    fexo.QARecords=FEMesh.QARecords;
end

fexo.InfoRecords=obj.rd_Info;
%%%%
[coord,nmap,names]=obj.rd_Nodes;
fexo.Nodes=FEMesh.Nodes(coord,nmap,names);
clear coord
%
fexo.CoordinateFrames=obj.rd_CoordFrames;

ids=obj.rd_BlockIDs;
fexo.Blocks=FEMesh.Blocks(ids);
for i=1:length(ids),
    [fexo.Blocks(i).Connectivity,fexo.Blocks(i).ElementType,fexo.Blocks(i).Name,fexo.Blocks(i).Status,fexo.Blocks(i).Attributes,fexo.Blocks(i).AttributesName]=obj.rd_Block(ids(i));
end

fexo.ElementMap=obj.rd_EMap;

nem=obj.rd_Nemesis;
if ~isempty(nem),
   fexo.Nemesis=nem;
end

%% read superelements if necessary
[kr,cr,mr,cb,otm,om,otme,oem,ddof,din,deig,dnc]=obj.rd_SuperElement;
if ~isempty(cb),
    fexo.SuperElement=FEMesh.SuperElement(kr,cr,mr,cb,otm,om,otme,oem,din,deig,dnc);
end
    
if flag<5,
    ids=obj.rd_NodesetIDs;
    fexo.Nodesets=FEMesh.Nodesets(ids);
    if ~isempty(ids),
        for i=1:length(ids),
            [fexo.Nodesets(i).Nodes,fexo.Nodesets(i).DistFactors,fexo.Nodesets(i).Name,fexo.Nodesets(i).Status]=obj.rd_Nodeset(ids(i));
        end
    end
else
    fexo.Nodesets=FEMesh.Nodesets([]);
end

if flag<4,
    ids=obj.rd_SidesetIDs;
    fexo.Sidesets=FEMesh.Sidesets(ids);
    if ~isempty(ids),
        for i=1:length(ids),
            [fexo.Sidesets(i).Elements,fexo.Sidesets(i).Sides,fexo.Sidesets(i).DistFactors,fexo.Sidesets(i).Status,fexo.Sidesets(i).Name]=obj.rd_Sideset(ids(i));
        end
    end
else
    fexo.Sidesets=FEMesh.Sidesets([]);
end

if flag<3 || flag == 6
    fexo.Time=obj.rd_Time;
else
    fexo.Time=[];
end

if flag<3,
    gnames=obj.rd_GlobalVarsNames;
    fexo.GlobalVars=FEMesh.GlobalVars(gnames);
    if ~isempty(gnames),
%         wh = waitbar(0,'Reading global variables');
%         waitInterval = 100;
        for i=1:length(gnames),
%             fexo.GlobalVars(i).Data=obj.rd_GlobalVar(gnames{i});
            fexo.GlobalVars(i).Data=obj.rd_GlobalVarByInd(i);
%             if ~mod(i,waitInterval)
%                 waitbar(i/length(gnames),wh,'Reading global variables');
%             end
        end
%         close(wh)
    end
else
    fexo.GlobalVars=FEMesh.GlobalVars([]);
end

if flag<2,
    nnames=obj.rd_NodalVarsNames;
    fexo.NodalVars=FEMesh.NodalVars(nnames);
    if ~isempty(nnames),
        for i=1:length(nnames),
            fexo.NodalVars(i).Data=obj.rd_NodalVar(nnames{i});
        end
    end
else
    fexo.NodalVars=FEMesh.NodalVars([]);
end


ids=fexo.Blocks.getBlockIDs;
if flag<1,
    enames=obj.rd_ElemVarsNames;
    fexo.ElemVars=FEMesh.ElemVars(ids,enames);
    if ~isempty(enames),
        for i=1:length(ids),
            for j=1:length(enames),
                fexo.ElemVars(i,j).Data=obj.rd_ElemVar(enames{j},ids(i));
            end
        end
    end
else
    fexo.ElemVars=FEMesh.ElemVars(ids,[]);
end

ids=fexo.Nodesets.getNodesetIDs;
if flag<1,
    nnames=obj.rd_NodesetVarsNames;
    fexo.NodesetVars=FEMesh.NodesetVars(ids,nnames);
    if ~isempty(nnames),
        for i=1:length(ids),
            for j=1:length(nnames),
                fexo.NodesetVars(i,j).Data=obj.rd_NodesetVar(nnames{j},ids(i));
            end
        end
    end
else
    fexo.NodesetVars=FEMesh.NodesetVars(ids,[]);
end

ids=fexo.Sidesets.getSidesetIDs;
if flag<1,
    nnames=obj.rd_SidesetVarsNames;
    fexo.SidesetVars=FEMesh.SidesetVars(ids,nnames);
    if ~isempty(nnames),
        for i=1:length(ids),
            for j=1:length(nnames),
                fexo.SidesetVars(i,j).Data=obj.rd_SidesetVar(nnames{j},ids(i));
            end
        end
    end
else
    fexo.SidesetVars=FEMesh.SidesetVars(ids,[]);
end
obj.close