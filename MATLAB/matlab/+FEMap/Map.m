classdef Map
    properties
        VariableTypeTo   %% Nodal or Element
        VariableTypeFrom
        ModelObjFrom
        ModelObjTo
    end
    methods
        function obj=Map(fexoFrom,fexoTo)
            if nargin>0,
                obj.ModelObjFrom=fexoFrom;
                obj.ModelObjTo=fexoTo;
            end
        end
        %%%%%%%%
        function obj=WriteNodal(obj,nvar,nname,time)
            obj.ModelObjTo=obj.ModelObjTo.AddNodalVar(nname,nvar,time);
        end
        function obj=WriteElem(obj,evar,ename,blkidx,nelem,ngauss,gp_id,time)
            jj=0;
            stride=0;
            
            for i=1:length(blkidx),
                %nelem=size(obj.ModelObjTo.Blocks(blkidx(i)).Connectivity,1);
                for j=1:ngauss(i),
                    obj=obj.AddElem(evar((jj+j):ngauss(i):(jj+nelem(i)*ngauss(i))), ...
                        ename,blkidx(i),gp_id(stride+j,:),ngauss(i),time);
                end
                stride=stride+ngauss(i);
                jj=jj+nelem(i)*ngauss(i);
            end
        end
        function obj=WriteElem2(obj,evar,ename,blkidx,nelem,ngauss,gp_id,time)
            %%%%%%
            ii=1;
            ng=[0;ngauss(:)];
            for i=1:length(blkidx),
                for j=1:ngauss(i),
                    idx=ii:ii+nelem(i)-1;
                    %evar=detJ*weights; % for volume
                    %%%%%%
                    ii=ii+nelem(i);
                    obj=obj.AddElem(evar(idx),ename,blkidx(i),gp_id(j+ng(i),:),ngauss(i),time);
                end
            end
        end
        function obj=AddElem(obj,evar,ename,blkidx,gaussid,ngauss,time)
            %% changed order of variables; ename and blkidx switched
            enameg=obj.mkEVarName(obj,blkidx,ename,gaussid,ngauss);
            
            obj.ModelObjTo=obj.ModelObjTo.AddElemVar(enameg,blkidx,evar,time);
        end
        function [ptxyz,evar]=getEVar(obj,blkidx,ename,tidx)
            %% Assumptions:  all similar elements eg hex20 have the same
            %% number of gauss points.
            %% TODO:  This code is terrible on memory.  Need to preallocate arrays
            
            ndim=size(obj.ModelObjFrom.Nodes.Coordinates,2);
            xc=obj.ModelObjFrom.Nodes.Coordinates(:,1);
            if ndim>1,
                yc=obj.ModelObjFrom.Nodes.Coordinates(:,2);
            else
                yc=zeros(size(xc));
            end
            if ndim>2,
                zc=obj.ModelObjFrom.Nodes.Coordinates(:,3);
            else
                zc=zeros(size(xc));
            end
            %%%
            ptxyz=[];%% Need to preallocate these arrays
            evar=[];
            for i=1:length(blkidx),
                names=obj.ModelObjFrom.ElemVars.getNames(blkidx(i));
                %idxe=strmatch(ename,names,'exact');
                idxe=find(strcmp(ename,names));
                
                strmatch=sprintf('%s_%s%d',ename, ...
                    upper(obj.ModelObjFrom.Blocks(blkidx(i)).ElementType(1:3)), ...
                    size(obj.ModelObjFrom.Blocks(blkidx(i)).Connectivity,2));
                %idxg=strmatch(sprintf('%s_%s%d',ename, ...
                %    upper(obj.ModelObjFrom.Blocks(blkidx(i)).ElementType(1:3)), ...
                %    size(obj.ModelObjFrom.Blocks(blkidx(i)).Connectivity,2)),names);
                idxg=find(strncmp(strmatch,names,length(strmatch)));
                
                ShpObj=ShapeFactory.CreateShape(obj.ModelObjFrom.Blocks(blkidx(i)).ElementType, ...
                    obj.ModelObjFrom.Blocks(blkidx(i)).Connectivity(1,:), ...
                    [1 1 1]);
                %%%%
                if ~isempty(idxg) %% Use Gauss data if it is available
                    %% get the gauss string
                    gcell=cell(length(idxg),1);
                    for j=1:length(idxg),
                        gcell{j}=names{idxg(j)}(end-ndim+1:end);
                    end
                    gstr=char(gcell);
                    ngpts=[0 0 0];
                    for j=1:ndim,
                        ngpts(j)=1+max(str2double(gstr(:,j)));
                    end
                    [gpts,wts,g_id]=ShpObj.GaussWeights(ngpts);
                    N=ShpObj.Interpolate(gpts);
                    ev=[];%%%  not the best approach but easy and quick
                    pxyz=[];  %% need to preallocate these arrays
                    for j=1:length(idxg),
                        idx=strmatch(sprintf('%d',g_id(j,:)),gcell);
                        ev=cat(1,ev,obj.ModelObjFrom.ElemVars(idxg(idx)).Data(:,tidx));
                        pxyz=cat(1,pxyz,[xc(obj.ModelObjFrom.Blocks(blkidx(i)).Connectivity)*N' ...
                            yc(obj.ModelObjFrom.Blocks(blkidx(i)).Connectivity)*N' ...
                            zc(obj.ModelObjFrom.Blocks(blkidx(i)).Connectivity)*N']);
                        %% here is where pxyz and ev are populated
                    end
                elseif ~isempty(idxe) %% use the centroid data if necessary
                    gpts=ShpObj.LocalCentroid;
                    N=ShpObj.Interpolate(gpts);
                    pxyz=[xc(obj.ModelObjFrom.Blocks(blkidx(i)).Connectivity)*N' ...
                        yc(obj.ModelObjFrom.Blocks(blkidx(i)).Connectivity)*N' ...
                        zc(obj.ModelObjFrom.Blocks(blkidx(i)).Connectivity)*N'];
                    ev=obj.ModelObjFrom.ElemVars(blkidx(i),idxe).Data(:,tidx);
                else   %% the variable wasn't defined in database
                    error('Variable %s not found in model',ename)
                end
                ptxyz=cat(1,ptxyz,pxyz);
                evar=cat(1,evar,ev);
            end
        end
        function [evar]=getEVarFrom(obj,ename,tidx)
            evar=[];
            for i=1:length(obj.ModelObjFrom.Blocks),
                names=obj.ModelObjFrom.ElemVars.getNames(i);
                %idx=strmatch(ename,names,'exact');
                idx=find(strcmp(ename,names));
                evar=cat(1,evar,obj.ModelObjFrom.ElemVars(i,idx).Data(:,tidx));
            end
        end
        
        %%%%
        function plot(obj,fromelem,toelem)
            h=obj.ModelObjFrom.plote(fromelem,[1 0 0]);
            obj.ModelObjTo.plote(toelem,[0 0 0],h);
            title(sprintf('%d From Elements(red), %d To Elements (black)',length(fromelem),length(toelem)))
        end
        function plotblk(obj,fromelem,toelem)
            a=238/255;
            b=138/255;
            [el,blk]=obj.ModelObjFrom.Blocks.LocalElemNum(fromelem);
            h=obj.ModelObjFrom.plot(blk,[0 a 0]);
          
            [el,blk]=obj.ModelObjTo.Blocks.LocalElemNum(toelem);
            %obj.ModelObjTo.plot(blk,[b b b],h);

            obj.ModelObjFrom.plote(fromelem,[1 0 0],h);
            obj.ModelObjTo.plote(toelem,[0 0 0],h);
            title(sprintf('%d "From" Elements(red), %d "To" Elements (black)',length(fromelem),length(toelem)))
            hold off
        end
        %%%%%%
    end
    methods (Abstract)
        %nvar=NodalVarMapFnc(obj,nodex,nodey,nodez)
        %evar=ElemVarMapFnc(obj,nodex,nodey,nodez)
    end
    methods (Static)
        function [enodexyz,nelem,gp_id,ng]=CalcElemCoords(ModelObj, ...
                blkidx,numgauss)
            %% ng is the number of gauss points per block
            %numgausspts is the number of final gauss point data you want
            %% in the exodus file
            %% numintpts is the number of points used to integrate over the
            %% volume in the subcell
            %%%
            nelem=zeros(length(blkidx),1);
            if ~numgauss,
                ngauss=zeros(length(blkidx),1);
            else
                ngauss=ones(length(blkidx),1);
            end
            for i=1:length(blkidx),
                nelem(i)=size(ModelObj.Blocks(blkidx(i)).Connectivity,1);
            end
            %%%%%
            dim=size(ModelObj.Nodes.Coordinates,2);
            enodexyz=zeros(sum(nelem.*ngauss),3);
           
            %%%%%
            gp_id=[];
            
            %%%
            ii=1;
            x=ModelObj.Nodes.Coordinates(:,1);
            y=ModelObj.Nodes.Coordinates(:,2);
            if dim==2,
                z=zeros(size(x));
            else
                z=ModelObj.Nodes.Coordinates(:,3);
            end
            ng=zeros(length(ngauss),1);
            for i=1:length(blkidx),
                ShpObj=ShapeFactory.CreateShape(ModelObj.Blocks(blkidx(i)).ElementType,ModelObj.Blocks(blkidx(i)).Connectivity,ModelObj.Nodes.Coordinates);
                
                if ngauss(i),
                    [gausspts,ju,wr_idx]=ShpObj.getIntPts(ShpObj.NumMassIntPts);
                    ng(i)=size(gausspts,1);
                else
                    [gausspts,ju,wr_idx]=ShpObj.getIntPts(1);
                    ng(i)=1;
                end
                gp_id=cat(1,gp_id,wr_idx);
                
                N=ShpObj.Interpolate(gausspts);
                %gausspts
                %x(ModelObj.Blocks(blkidx(i)).Connectivity)*N'
                
                enodexyz(ii:ii+nelem(i)*ng(i)-1,:)=[reshape(x(ModelObj.Blocks(blkidx(i)).Connectivity)*N',nelem(i)*ng(i),1) ...
                    reshape(y(ModelObj.Blocks(blkidx(i)).Connectivity)*N',nelem(i)*ng(i),1) ...
                    reshape(z(ModelObj.Blocks(blkidx(i)).Connectivity)*N',nelem(i)*ng(i),1)];
                
                ii=ii+nelem(i)*ng(i);
                
            end
        end
        function enameg=mkEVarName(obj,blkidx,ename,gaussid,ngauss)
            if ngauss==1,
                enameg=ename;
            else
                elemtype=obj.ModelObjTo.Blocks(blkidx).ElementType;
                nnodes=size(obj.ModelObjTo.Blocks(blkidx).Connectivity,2);
                switch lower(elemtype(1:3)),
                    case 'hex',
                        etype='HEX';
                    case 'qua',
                        etype='QUA';
                    case 'she',
                        etype='SHE';
                    case 'tet',
                        etype='TET';
                    case 'wed'
                        etype='WED';
                    otherwise
                        error('Unknown Element Type %s:Map',elemtype)
                end
                switch size(gaussid,2),
                    case 1
                        enameg=sprintf('%s_%s%d_GP%d',ename,etype,nnodes,gaussid(1));
                    case 2
                        enameg=sprintf('%s_%s%d_GP%d%d',ename,etype,nnodes,gaussid(1),gaussid(2));
                    case 3
                        enameg=sprintf('%s_%s%d_GP%d%d%d',ename,etype,nnodes,gaussid(1),gaussid(2),gaussid(3));
                end
            end
        end
        function nodexyz=ElemCoord(ModelObj,blkidx,ShpObj,intpt)
            x=ModelObj.Nodes.Coordinates(:,1);
            y=ModelObj.Nodes.Coordinates(:,2);
            
            switch size(ModelObj.Nodes.Coordinates,2)
                case 2
                    z=zeros(size(x));
                case 3
                    z=ModelObj.Nodes.Coordinates(:,3);
                otherwise
                    error('unknown dimension, %d',size(ModelObj.Nodes.Coordinates,2))
            end
            %%
            switch length(intpt),
                case 2
                    N=ShpObj.Interpolate(intpt(1),intpt(2),0);
                case 3
                    N=ShpObj.Interpolate(intpt(1),intpt(2),intpt(3));
                otherwise
                    error('Weird Number of gauss coordinates')
            end
            
            switch size(ModelObj.Nodes.Coordinates,2),
                case 2
                    nodexyz=[x(ModelObj.Blocks(blkidx).Connectivity)*N', ...
                        y(ModelObj.Blocks(blkidx).Connectivity)*N', ...
                        zeros(size(ModelObj.Blocks(blkidx).Connectivity,1),1)];
                    %%%%%%
                case 3
                    nodexyz=[x(ModelObj.Blocks(blkidx).Connectivity)*N', ...
                        y(ModelObj.Blocks(blkidx).Connectivity)*N', ...
                        z(ModelObj.Blocks(blkidx).Connectivity)*N'];
                    %%%%%%
                otherwise
                    error('Unknown Dimension')
            end
        end
        function [int_evar,tot_vol,cint_evar,ctot_vol]=IntBlock(ModelObj,blkidx,ename,tidx)
            %% Assumptions:  all similar elements eg hex20 have the same
            %% number of gauss points.
            %%
            %%  cint_evar and ctot_vol are the corrected volumes for data
            %%  from ITS.  The code corrects the integration based upon the
            %%  variable VOL which gives the true volume of material in
            %%  each element.
            %%
            %% TODO:  This code is terrible on memory.  Need to preallocate arrays
            if length(blkidx)>1,
                error('IntBlock only operates on one block at a time')
            end
            ndim=size(ModelObj.Nodes.Coordinates,2);
            
            %%%
            
            names=ModelObj.ElemVars.getNames;
            
            %idxe=strmatch(ename,names,'exact');
            idxe=find(strcmp(ename,names));
            strmatch=sprintf('%s_%s%d',ename, ...
                upper(ModelObj.Blocks(blkidx).ElementType(1:3)), ...
                size(ModelObj.Blocks(blkidx).Connectivity,2));
            %idxg=strmatch(sprintf('%s_%s%d',ename, ...
            %    upper(ModelObj.Blocks(blkidx).ElementType(1:3)), ...
            %    size(ModelObj.Blocks(blkidx).Connectivity,2)),names);
            idxg=find(strncmp(strmatch,names,length(strmatch)));

            ShpObj=ShapeFactory.CreateShape(ModelObj.Blocks(blkidx).ElementType, ...
                ModelObj.Blocks(blkidx).Connectivity, ...
                ModelObj.Nodes.Coordinates);
            %%%%
            if ~isempty(idxg) %% Use Gauss data if it is available
                %% get the gauss string
                gcell=cell(length(idxg),1);
                for j=1:length(idxg),
                    gcell{j}=names{idxg(j)}(end-ndim+1:end);
                end
                gstr=char(gcell);
                ngpts=[0 0 0];
                for j=1:ndim,
                    ngpts(j)=1+max(str2double(gstr(:,j)));
                end
                %[gpts,wts,g_id]=ShpObj.GaussWeights(ngpts);
                [gpts,wts,g_id]=ShpObj.getIntPts(ShpObj.NumMassIntPts);
                   
                nelem=size(ModelObj.Blocks(blkidx).Connectivity,1);
                int_evar=0;
                tot_vol=0;
                for j=1:length(idxg),
                    strmatch=sprintf('%d',g_id(j,:));
                    %idx=strmatch(sprintf('%d',g_id(j,:)),gcell);
                    idx=find(strncmp(strmatch,gcell,length(strmatch)));
                    if isempty(idx),
                        error('Something is wrong')
                    end
                    [detJ]=ShpObj.calcDetJ(repmat(gpts(j,:),nelem,1),1:nelem);
                    int_evar=int_evar+sum(wts(j)*ModelObj.ElemVars(blkidx,idxg(idx)).Data(:,tidx).*detJ);
                    tot_vol=tot_vol+sum(wts(j)*detJ);
                end
        
            elseif ~isempty(idxe) %% use the centroid data if necessary
                vol=abs(ShpObj.Volume);
                int_evar=sum(ModelObj.ElemVars(blkidx,idxe).Data(:,tidx).*vol);
                tot_vol=sum(vol);
            else   %% the variable wasn't defined in database
                error('Variable %s not found in model',ename)
            end
            if ModelObj.ElemVars.isElemVar('Vol')
                %idxv=strmatch('VOL',names,'exact');
                %idxv=strmatch('VOL',names);
                idxv=find(strncmp('VOL',names,length('VOL')));
                volc=ModelObj.ElemVars(blkidx,idxv).Data(:,tidx);
                cint_evar=sum(ModelObj.ElemVars(blkidx,idxe).Data(:,tidx).*volc);
                ctot_vol=sum(volc);
            else
                cint_evar=[];
                ctot_vol=[];
            end
        end
        function [evar,coord]=GaussData(ModelObj,blkidx,ename,tidx)
            
            %% TODO:  This code is terrible on memory.  Need to preallocate arrays
            if length(blkidx)>1,
                error('GaussData only operates on one block at a time')
            end
            ndim=size(ModelObj.Nodes.Coordinates,2);
            
            %%%
            
            names=ModelObj.ElemVars.getNames(blkidx);
            %idxe=strmatch(ename,names,'exact');
            idxe=find(strcmp(ename,names));
            strmatch=sprintf('%s_%s%d',ename, ...
                upper(ModelObj.Blocks(blkidx).ElementType(1:3)), ...
                size(ModelObj.Blocks(blkidx).Connectivity,2));
            %idxg=strmatch(sprintf('%s_%s%d',ename, ...
            %    upper(ModelObj.Blocks(blkidx).ElementType(1:3)), ...
            %    size(ModelObj.Blocks(blkidx).Connectivity,2)),names);
            idxg=find(strncmp(strmatch,names,length(strmatch)));
            ShpObj=ShapeFactory.CreateShape(ModelObj.Blocks(blkidx).ElementType, ...
                ModelObj.Blocks(blkidx).Connectivity, ...
                ModelObj.Nodes.Coordinates);
            %%%%
            if ~isempty(idxg) %% Use Gauss data if it is available
                %% get the gauss string
                gcell=cell(length(idxg),1);
                for j=1:length(idxg),
                    gcell{j}=names{idxg(j)}(end-ndim+1:end);
                end
                gstr=char(gcell);
                ngpts=[0 0 0];
                for j=1:ndim,
                    ngpts(j)=1+max(str2double(gstr(:,j)));
                end
                %[gpts,wts,g_id]=ShpObj.GaussWeights(ngpts);
                [gpts,wts,g_id]=ShpObj.getIntPts(ShpObj.NumMassIntPts);
                   
                nelem=size(ModelObj.Blocks(blkidx).Connectivity,1);
                coord=[];
                evar=[];
                for j=1:length(idxg),
                    idx=strmatch(sprintf('%d',g_id(j,:)),gcell);
                    coord=cat(1,coord,ShpObj.Local2Global(repmat(gpts(j,:),nelem,1),1:nelem));
                    evar=cat(1,evar,ModelObj.ElemVars(blkidx,idxg(idx)).Data(:,tidx));
                end
        
            elseif ~isempty(idxe) %% use the centroid data if necessary
                
                evar=ModelObj.ElemVars(blkidx,idxe).Data(:,tidx);
                coord=ShpObj.GlobalCentroid;
            else   %% the variable wasn't defined in database
                error('Variable %s not found in model',ename)
            end
            
        end
        %%
        function fexoTo=CorrEvar(fexoTo,fexoFrom,ename,blkidxTo,blkidxFrom,tidxTo,tidxFrom)
            % scales the fexoTo variable evar such that the volume ingrated 
            %   evar is equal in both (energy is typically used)
            %
            
            for i=1:length(tidxTo),
                tot_from_evar=0;
                for j=1:length(blkidxFrom), % figure out total energy on original mesh
                    [int_evar,tv]=IntBlock(fexoFrom,blkidx,ename,fexoFrom.Time(tidxFrom(i)));
                    tot_from_evar=tot_from_evar+int_evar*tv;
                end
                tot_to_evar=0;
                for j=1:length(blkidxTo), % figure out total energy on original mesh
                    [int_evar,tv]=IntBlock(fexoTo,blkidx,ename,fexoTo.Time(tidxTo(i)));
                    tot_to_evar=tot_to_evar+int_evar*tv;
                end
                
                g=(tot_from_evar/tot_to_evar);
                for j=1:length(blkidxTo),
                    fexoTo.ElemVars(i).Data(:,tidxTo(i))=g*fexoTo.ElemVars(i).Data(:,tidxTo(i));
                end
            end
        end
        function fexo=CorrEnergy(fexo,ename,tidx)
            evarFrom=fexo.getElemVar(ename,tidx);

            time=fexo.Time(tidx);
            %%
            
            [ITSvolFrom]=fexo.getElemVar('VOL',tidx);
            volFrom=fexo.Volume_elem;
            evarFrom=(ITSvolFrom./volFrom).*evarFrom;
            
            jj=0;
            for i=1:length(fexo.Blocks),
                nelem=size(fexo.Blocks(i).Connectivity,1);
                
                fexo=fexo.AddElemVar(ename,i,evarFrom((jj+1):(jj+nelem)),time);
                jj=jj+nelem;
            end
            
        end
    end
end