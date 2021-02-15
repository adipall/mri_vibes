classdef NodalInterp < FEMap.Map
    properties
       Points
       PointsNotMapped
       LocalCoord
       GaussIDTo
       NumGaussPtsTo
       Elements
       Blocks
       NumberElementsperBlock
       BlkIdxTo
       BlkIdxFrom
    end
    
    methods
        function obj=NodalInterp(fexoFrom,fexoTo)
            if nargin>0,
                obj.ModelObjTo=fexoTo;
                obj.ModelObjFrom=fexoFrom;
            end
        end
        function obj=N2EMap(obj,blkidxfrom,blkidxto,ngauss,varnamefrom,varnameto,tidx)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%  Start heere look at ngauss
            obj.BlkIdxFrom=blkidxfrom;
            obj.BlkIdxTo=blkidxto;
            [enodexyz,obj.NumberElementsperBlock,obj.GaussIDTo,ng]=obj.CalcElemCoords(obj.ModelObjTo, ...
                obj.BlkIdxTo,ngauss);
            %%%
            obj.NumGaussPtsTo=ng;
            %%%%
 
            SrchObj=FESearch.ElementNode(obj.ModelObjFrom,obj.BlkIdxFrom);
            %%%
            [obj.LocalCoord,obj.Elements,obj.Blocks, ...
                obj.Points,obj.PointsNotMapped]=SrchObj.LocalCoord(enodexyz);
            %%%
            obj.VariableTypeTo='Element';
            obj.VariableTypeFrom='Node';
            for i=1:length(varnamefrom),
                if iscell(varnamefrom),
                    vnamefrom=varnamefrom{i};
                    vnameto=varnameto{i};
                else
                    vnamefrom=varnamefrom;
                    vnameto=varnameto;
                end
                obj=obj.Interp(vnamefrom,vnameto,tidx);
            end
        end
        function obj=N2NMap(obj,blkidxfrom,blkidxto,varnamefrom,varnameto,tidx)
            %if nargin<2,
            %    blkidxto=1:length(obj.ModelObjTo.Blocks);
            %end
            %if nargin<3,
            %    blkidxfrom=1:length(obj.ModelObjFrom.Blocks);
            %end
            obj.BlkIdxFrom=blkidxfrom;
            obj.BlkIdxTo=blkidxto;
            nodes=obj.ModelObjTo.Blocks.NodesInBlocks(blkidxto);
            
            SrchObj=FESearch.ElementNode(obj.ModelObjFrom,blkidxfrom);
            [obj.LocalCoord,obj.Elements,obj.Blocks, ...
                obj.Points,obj.PointsNotMapped]=SrchObj.LocalCoord(obj.ModelObjTo.Nodes.Coordinates(nodes,:));
            obj.Points=nodes(obj.Points);
            obj.PointsNotMapped=nodes(obj.PointsNotMapped);
            %%%%%
            obj.VariableTypeTo='Node';
            obj.VariableTypeFrom='Node';
            
            if iscell(varnamefrom),
                for i=1:length(varnamefrom),
                    vnamefrom=varnamefrom{i};
                    vnameto=varnameto{i};
                    obj=obj.Interp(vnamefrom,vnameto,tidx);
                end
            else
                vnamefrom=varnamefrom;
                vnameto=varnameto;
                obj=obj.Interp(vnamefrom,vnameto,tidx);
            end
        end
        %%%%%%%%%%%
        function obj=Interp(obj,varnamefrom,varnameto,tidx)
            nvar=obj.ModelObjFrom.NodalVars.getnvar(varnamefrom);

            var=zeros(length(unique([obj.Points(:);obj.PointsNotMapped(:)])),1);

            for i=1:length(obj.BlkIdxFrom),
                blk=obj.BlkIdxFrom(i);
                idx=find(i==obj.Blocks);
                ShpObj=ShapeFactory.CreateShape(obj.ModelObjFrom.Blocks(blk).ElementType, ...
                    obj.ModelObjFrom.Blocks(blk).Connectivity, ...
                    obj.ModelObjFrom.Nodes.Coordinates);
                N=ShpObj.Interpolate(obj.LocalCoord(idx,:));
                var(obj.Points(idx))=sum(N.*nvar(obj.ModelObjFrom.Blocks(blk).Connectivity(obj.Elements(idx),:)),2); 
            end
            %%%%%%
            switch obj.VariableTypeTo
                case 'Node'
                    obj=WriteNodal(obj,var,varnameto,obj.ModelObjFrom.Time(tidx));
                case 'Element'
                    obj=obj.WriteElem(var,varnameto,obj.BlkIdxTo,obj.NumberElementsperBlock, ...
                        obj.NumGaussPtsTo,obj.GaussIDTo,obj.ModelObjFrom.Time(tidx));
                
            end
        end
    end
end