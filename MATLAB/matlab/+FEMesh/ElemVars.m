classdef ElemVars < handle
    properties
        BlockID
        Name
        Data
    end
    
    methods
        function obj=ElemVars(blkids,names)
            if nargin>0,
                if nargin>1 && ~isempty(names),
                    nv=length(names);
                else
                    nv=1;
                end
                for i=1:length(blkids),
                    for j=1:nv,
                        obj(i,j)=FEMesh.ElemVars;
                    end
                end
                for i=1:length(blkids),
                    for j=1:nv,
                        obj(i,j).BlockID=blkids(i);
                        if nargin>1 && ~isempty(names),
                            obj(i,j).Name=names{j};
                        end
                    end
                end

            end
        end
        function display(obj)
            [n,m]=size(obj);
            fprintf(' BlockID      Variable Name\n')
            for i=1:n,
                for j=1:m,
                    if j==1,
                        str=sprintf('   %d         ',obj(i,j).BlockID);
                    end
                    if ~isempty(obj(i,j).Data)
                        str=[str,sprintf('%s',obj(i,j).Name)];
                        if j~=m,
                            str=[str,', '];
                        end
                    end
                end
                disp(str)
            end
        end
        
        function names=getNames(obj,blkidx)
            if nargin==2,
                ii=1;
                
                for i=1:size(obj,2),
                    if ~isempty(obj(blkidx,i).Data),
                        names{ii}=obj(blkidx,i).Name;
                        ii=ii+1;
                    end
                end
                if ii==1,
                    error('No element variables found for block index %d',blkidx)
                end
            else
                names=cell(size(obj,2),1);
                for i=1:size(obj,2),
                    names{i}=obj(1,i).Name;
                end
            end
        end
        function blkids=getBlockIDs(obj)
            blkids=zeros(size(obj,1),1);
            for i=1:size(obj,1),
                blkids(i)=obj(i,1).BlockID;
            end
        end
        function [idx]=id2idx(obj,id)
            ids=[obj(:,1).BlockID];
            idx=zeros(size(id));
            for i=1:length(id),
                vi=find(ids==id(i));
                if ~isempty(vi),
                    idx(i)=vi;
                else
                    error('Block ID %d not found',id(i))
                end
            end
        end
        function flag=isElemVar(obj,ename)
            flag=0;
            for j=1:size(obj,2),
                if strcmpi(obj(1,j).Name,ename)
                    flag=j;
                    break;
                end
            end
        end
    end
end
