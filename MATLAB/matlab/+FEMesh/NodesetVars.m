classdef NodesetVars < handle
    properties
        NodesetID
        Name
        Data
    end
    
    methods
        function obj=NodesetVars(nsids,names)
            if nargin>0,
                if nargin>1 && ~isempty(names),
                    nv=length(names);
                else
                    nv=1;
                end
                for i=1:length(nsids),
                    for j=1:nv,
                        obj(i,j)=FEMesh.NodesetVars;
                    end
                end
                for i=1:length(nsids),
                    for j=1:nv,
                        obj(i,j).NodesetID=nsids(i);
                        if nargin>1 && ~isempty(names),
                            obj(i,j).Name=names{j};
                        end
                    end
                end

            end
        end
        function display(obj)
            [n,m]=size(obj);
            fsprintf(' NodesetID      Variable Name\n')
            for i=1:n,
                for j=1:m,
                    if j==1,
                        str=sprintf('   %d         ',obj(i,j).NodesetID);
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
        
        function names=getNames(obj,nsidx)
            if nargin==2,
                ii=1;
                
                for i=1:size(obj,2),
                    if ~isempty(obj(nsidx,i).Data),
                        names{ii}=obj(nsidx,i).Name;
                        ii=ii+1;
                    end
                end
                if ii==1,
                    error('No element variables found for nodeset index %d',nsidx)
                end
            else
                names=cell(size(obj,2),1);
                for i=1:size(obj,2),
                    names{i}=obj(1,i).Name;
                end
            end
        end
        function nsids=getNodesetIDs(obj)
            nsids=zeros(size(obj,1),1);
            for i=1:size(obj,1),
                nsids(i)=obj(i,1).NodesetID;
            end
        end
        function flag=isNodesetVar(obj,nname)
            flag=0;
            for j=1:size(obj,2),
                if strcmpi(obj(1,j).Name,nname)
                    flag=j;
                    break;
                end
            end
        end
    end
end
