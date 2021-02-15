classdef ElemMap < FEMap.Map
    properties
        AbsoluteSearchTolerance = []
        LocalCoordinateTolerance = []
        TransferModelFrom
        TransferModelTo
        BlkIdxTo
        NumGaussPtsTo
        GaussIDTo
        NumElementsTo
    end
    methods
        function obj=ElemMap(fexoFrom,fexoTo)
            if nargin>0,
                obj.ModelObjTo=fexoTo;
                obj.ModelObjFrom=fexoFrom; 
            end
        end
        function obj=setBlockIdx(obj,blkidxfrom,blkidxto,gauss)
            %  if gauss==1 map to gauss points
            %   if gauss==0 map to centroid
            %
            obj.TransferModelFrom=[];
            obj.TransferModelFrom=obj.ModelObjFrom.GaussMeshRefinement(blkidxfrom, ...
                ones(length(blkidxfrom),1),'gauss');
            obj.TransferModelFrom.ElemVars=obj.ModelObjFrom.ElemVars(blkidxfrom,:);
            obj.TransferModelFrom.Time=obj.ModelObjFrom.Time;
            nelem=0;
            gid=[];
            ngaussto=ones(length(blkidxto),1);
            obj.NumGaussPtsTo=zeros(length(blkidxto),1);
            if gauss==0,
                obj.TransferModelTo=obj.ModelObjTo.Copy.ModelReduce(blkidxto);
                gid=zeros(length(blkidxto),3);
                obj.NumGaussPtsTo=ones(length(blkidxto),1);
            else
                for i=1:length(blkidxto),
                    Sobj=ShapeFactory.CreateShape(obj.ModelObjTo.Blocks(blkidxto(i)).ElementType, ...
                        obj.ModelObjTo.Blocks(blkidxto(i)).Connectivity(i,:),[]);
                    if gauss==1,
                        [ju1,ju2,gid1]=Sobj.getIntPts(Sobj.NumMassIntPts);
                        ngaussto(i)=Sobj.NumMassIntPts;
                    else
                        gid1=[0 0 0];
                    end
                    gid=cat(1,gid,gid1);
                    obj.NumGaussPtsTo(i)=size(gid1,1);
                end
                obj.TransferModelTo=[];
                obj.TransferModelTo=obj.ModelObjTo.GaussMeshRefinement(blkidxto, ...
                    ngaussto,'gauss');
            end
            for i=1:length(blkidxto),
                nelem=nelem+size(obj.TransferModelTo.Blocks(i).Connectivity,1);
            end
            obj.GaussIDTo=gid;
            obj.NumElementsTo=nelem;
            obj.BlkIdxTo=blkidxto;
        end
        function [obj,Ielem,vol]=E2EMap(obj,ename,tidx,type,scale)
            if nargin<5,
                scale=1;
            end
            % use scale to better condition the problem if necessary
            switch type
                case 1 %% 
                    %% Rashid, International Journal of Numerical Methods in
                    %% Engineering, vol 55, pp 431-450, 2002.
                    blkidxfrom=1:length(obj.TransferModelFrom.Blocks);
                    blkidxto=1:length(obj.TransferModelTo.Blocks);
                    Sobj=FESearch.ElementElement(obj.TransferModelTo,blkidxto,scale);
                   
                    [vol,Ielem,truevol2]=Sobj.ElementIntersectionRashid(obj.TransferModelFrom,blkidxfrom);
                    
                    obj=obj.VolWeightEvars1(ename,vol,truevol2,Ielem,tidx);
                case 2 %%
                    %% Rashid, International Journal of Numerical Methods in
                    %% Engineering, vol 55, pp 431-450, 2002.
                    blkidxfrom=1:length(obj.TransferModelFrom.Blocks);
                    blkidxto=1:length(obj.TransferModelTo.Blocks);
                    Sobj=FESearch.ElementElement(obj.TransferModelFrom,blkidxfrom);
                    [vol,Ielem]=Sobj.ElementIntersectionRashid(obj.TransferModelTo,blkidxto);
                    
                    obj=obj.VolWeightEvars2(ename,vol,Ielem,tidx);
                case 3 %%
                    %% Rashid, International Journal of Numerical Methods in
                    %% Engineering, vol 55, pp 431-450, 2002.
                    blkidxfrom=1:length(obj.TransferModelFrom.Blocks);
                    blkidxto=1:length(obj.TransferModelTo.Blocks);
                    Sobj=FESearch.ElementElement(obj.TransferModelTo,blkidxto);
                    
                    [vol,Ielem,truevol2,volFrom]=Sobj.ElementIntersectionRashid(obj.TransferModelFrom,blkidxfrom);
                    %[vol,Ielem,truevol2,volFrom]=Sobj.ElementIntersectionClark(obj.TransferModelFrom,blkidxfrom);
                    obj=obj.VolWeightEvars3(ename,vol,truevol2,volFrom,Ielem,tidx);
                case 4 %%
                   %% Rashid, International Journal of Numerical Methods in
                    %% Engineering, vol 55, pp 431-450, 2002.
                    blkidxfrom=1:length(obj.TransferModelFrom.Blocks);
                    blkidxto=1:length(obj.TransferModelTo.Blocks);
                    Sobj=FESearch.ElementElement(obj.TransferModelFrom,blkidxfrom);
                    
                    %[vol,Ielem,truevolfrom,volto]=Sobj.ElementIntersectionRashid(obj.TransferModelTo,blkidxto,scale);
                    [vol,Ielem,truevolfrom,volto]=Sobj.ElementIntersectionClark(obj.TransferModelTo,blkidxto,scale);
                    obj=obj.VolWeightEvars4(ename,vol,truevolfrom,volto,Ielem,tidx);
            end
            obj.VariableTypeTo='Element';
            obj.VariableTypeFrom='Element';
        end
        function obj=WriteElem(obj,evar,ename,time)
            jj=0;
            stride=0;
            obj.NumGaussPtsTo
            for i=1:length(obj.BlkIdxTo),
                nelem=size(obj.ModelObjTo.Blocks(obj.BlkIdxTo(i)).Connectivity,1);
                for j=1:obj.NumGaussPtsTo(i),  
                    obj=obj.AddElem(evar((jj+j):obj.NumGaussPtsTo(i):(jj+nelem*obj.NumGaussPtsTo(i))), ...
                        ename,obj.BlkIdxTo(i),obj.GaussIDTo(stride+j,:),obj.NumGaussPtsTo(i),time);
                end
                stride=stride+obj.NumGaussPtsTo(i);
                jj=jj+nelem*obj.NumGaussPtsTo(i);
            end
        end
        function [evar]=getEVar(obj,ename,tidx)
            evar=[];
            for i=1:length(obj.TransferModelFrom.Blocks),
                names=obj.TransferModelFrom.ElemVars.getNames(i);
                idx=strmatch(ename,names,'exact');
                if ~isempty(idx),
                    evar=cat(1,evar,obj.TransferModelFrom.ElemVars(i,idx).Data(:,tidx));
                else
                    error('Variable %s not found',ename)
                end
            end
        end
            
    end
    methods (Access=private)
        function obj=VolWeightEvars1(obj,ename,vol,truevol2,Ielem,tidx)
            %%

            [evarFrom]=obj.getEVar(ename,tidx);
            
            evar=zeros(sum(obj.NumElementsTo),1);
            volest=zeros(sum(obj.NumElementsTo),1);

            for i=1:length(Ielem),
                if ~isempty(Ielem{i})
                    evar(Ielem{i})=evar(Ielem{i})+vol{i}*evarFrom(i);
                    volest(Ielem{i})=volest(Ielem{i})+vol{i};
                end
                %vol1est(i)=sum(vol{i});
            end
            evar=evar./truevol2;
            err=(volest-truevol2)./truevol2;
            %[volest truevol2 err]
            %[sum(volest) sum(truevol2)]
            %max(abs(vol1est-truevol1)./truevol1)
            obj=WriteElem(obj,evar,ename,obj.ModelObjFrom.Time(tidx));
            %
            obj=WriteElem(obj,err,'volerr',obj.ModelObjFrom.Time(tidx));
        end
        function obj=VolWeightEvars2(obj,ename,vol,Ielem,tidx)
            %%
            [evarFrom]=obj.getEVar(ename,tidx);
            evar=zeros(sum(obj.NumElementsTo),1);
            volest=zeros(sum(obj.NumElementsTo),1);
            %vol1est=zeros(length(Ielem),1);
            for i=1:length(Ielem),
                if ~isempty(Ielem{i}),
                    evar(i)=sum(vol{i}.*evarFrom(Ielem{i}));
                    volest(i)=sum(vol{i});
                    %vol1est(i)=sum(vol{i});
                end
            end
            
            evar=evar./volest;
            obj=WriteElem(obj,evar,ename,obj.ModelObjFrom.Time(tidx));
            %max(abs(vol1est-truevol1)./truevol1)
        end
        function obj=VolWeightEvars3(obj,ename,vol,truevol2,volFrom,Ielem,tidx)
            %%  Same as 1 except uses the var VOL from the exodus to
            %%  calculate the load.  This is added by ITS to account for
            %%  nonconformal elements
            [evarFrom]=obj.getEVar(ename,tidx);
            %%
            if ~obj.TransferModelFrom.ElemVars.isElemVar('VOL'),
                error('VOL Variable from ITS is not present')
            end
            [ITSvolFrom]=obj.getEVar('VOL',tidx);
            
            evarFrom=(ITSvolFrom./volFrom).*evarFrom;
            
            %%
            evar=zeros(sum(obj.NumElementsTo),1);
            volest=zeros(sum(obj.NumElementsTo),1);
            %size(evar)
            for i=1:length(Ielem),
                if ~isempty(Ielem{i})
                    %Ielem{i}
                    evar(Ielem{i})=evar(Ielem{i})+vol{i}*evarFrom(i);
                    volest(Ielem{i})=volest(Ielem{i})+vol{i};
                end
                %vol1est(i)=sum(vol{i});
            end
            evar=evar./truevol2;
            err=(volest-truevol2)./truevol2;
            %[volest truevol2 err]
            %[sum(volest) sum(truevol2)]
            %max(abs(vol1est-truevol1)./truevol1)
            obj=WriteElem(obj,evar,ename,obj.ModelObjFrom.Time(tidx));
            %
            %obj=WriteElem(obj,err,'volerr',obj.ModelObjFrom.Time(tidx));
           
        end
        function obj=VolWeightEvars4(obj,ename,vol,truevolfrom,volto,Ielem,tidx)
            %%  Same as 2 except uses the var VOL from the exodus to
            %%  calculate the load.  This is added by ITS to account for
            %%  nonconformal elements
            %%
            [evarFrom]=obj.getEVar_simp(ename,tidx);
            %%
            if ~obj.TransferModelFrom.ElemVars.isElemVar('VOL'),
                error('VOL Variable from ITS is not present')
            end
            [volFrom]=obj.getEVar_simp('VOL',tidx);
            
            evarFrom=(volFrom./truevolfrom).*evarFrom;
            %%
            evar=zeros(sum(obj.NumElementsTo),1);
            volest=zeros(sum(obj.NumElementsTo),1);
            %vol1est=zeros(length(Ielem),1);
            for i=1:length(Ielem),
                if ~isempty(Ielem{i}),
                    evar(i)=sum(vol{i}.*evarFrom(Ielem{i}));
                    volest(i)=sum(vol{i});
                    %vol1est(i)=sum(vol{i});
                end
            end
            err=(volest-volto)./volto;
            evar=evar./volto;
            obj=WriteElem(obj,evar,ename,obj.ModelObjFrom.Time(tidx));
            obj=WriteElem(obj,err,'volerr',obj.ModelObjFrom.Time(tidx));
            %max(abs(vol1est-truevol1)./truevol1)
        end
    end
    
end