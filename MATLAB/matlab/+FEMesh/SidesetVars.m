classdef SidesetVars < handle
    properties
        SidesetID
        Name
        Data
    end
    
    methods
        function obj=SidesetVars(ssids,names)
            if nargin>0,
                if nargin>1 && ~isempty(names),
                    nv=length(names);
                else
                    nv=1;
                end
                for i=1:length(ssids),
                    for j=1:nv,
                        obj(i,j)=FEMesh.SidesetVars;
                    end
                end
                for i=1:length(ssids),
                    for j=1:nv,
                        obj(i,j).SidesetID=ssids(i);
                        if nargin>1 && ~isempty(names),
                            obj(i,j).Name=names{j};
                        end
                    end
                end

            end
        end
        function display(obj)
            [n,m]=size(obj);
            fprintf(sprintf(' SidesetID      Variable Name\n'))
            for i=1:n,
                for j=1:m,
                    if j==1,
                        str=sprintf('   %d         ',obj(i,j).SidesetID);
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
                    error('No variables found for sideset index %d',nsidx)
                end
            else
                names=cell(size(obj,2),1);
                for i=1:size(obj,2),
                    names{i}=obj(1,i).Name;
                end
            end
        end
        function nsids=getSidesetIDs(obj)
            nsids=zeros(size(obj,1),1);
            for i=1:size(obj,1),
                ssids(i)=obj(i,1).SidesetID;
            end
        end
        function flag=isSidesetVar(obj,nname)
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
