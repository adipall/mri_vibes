classdef Blocks < handle
    properties
        ID
        ElementType
        Connectivity
        Name
        Attributes %%% nelems x nattributes
        AttributesName
        Status  = 1
    end
    properties (Dependent)
        NumElementsBlk
    end
    methods
        function obj=Blocks(ids)
            if nargin>0,
                n=length(ids);
                %% initializing obj arrays
                for i=1:n,
                    obj(i).ID=ids(i);
                end
            end
        end
        function display(obj)
            disp(sprintf('Index     ID       Name       NumElements      ElementType     NumNodesPerElement'))
            for i=1:length(obj),
                disp(sprintf('%d         %d       %s            %d                %s              %d',i,obj(i).ID, ...
                deblank(obj(i).Name'),obj(i).NumElementsBlk,obj(i).ElementType,size(obj(i).Connectivity,2)))
            end 
        end
        function set.Connectivity(obj,conn)
            obj.Connectivity=conn;
        end
        function num=numelem(obj)
            num=0;
            for i=1:length(obj),
                num=num+size(obj(i).Connectivity,1);
            end
        end
        function elemnum=GlobalElemNum(obj,blkidx,elem)
            start=[0 cumsum([obj.NumElementsBlk])];
            nelem=obj(blkidx).numelem;
            if elem>obj(blkidx).numelem,
                error('Maximum number of elements in Block %d is %d\n',blkidx,nelem);
            end
            elemnum=start(blkidx)+elem;
        end
        function [elemnum,blkidx]=LocalElemNum(obj,elem)
            start=[0 cumsum([obj.NumElementsBlk])];
            elemnum=zeros(length(elem),1);
            blkidx=zeros(length(elem),1);
            for i=1:length(elem)
                I=sort([start(:);elem(i)]);
                idx=find(I==elem(i),1,'first');
                elemnum(i)=elem(i)-start(idx-1);
                blkidx(i)=idx-1;
            end
        end
        function [node_elem,nodes]=inv_connect(obj)
            %
            % Determine how many nodes there are
            maxnode=-1;
            minnode=9e16;
            for i=1:length(obj),
                maxnode=max(max(max(obj(i).Connectivity)),maxnode);
                minnode=min(min(min(obj(i).Connectivity)),minnode);
            end
            node_elem=cell(maxnode-minnode+1,1);
            for i=1:length(obj),
                conn=obj(i).Connectivity;
                num_nodes_per_element=size(conn,2);
                %
                conn=conn';
                [scon,idx]=sort(conn(:));
                
                [nodes,idx2]=unique(scon,'last');
                %% convert to global indices on the elements
                %%idx2=obj.GlobalElemNum(i,idx2);
                idx=idx-1;
                
                idx=floor(idx/num_nodes_per_element)+1; % these are the elements corresponding to the nodes in scon
                idx(~idx)=1;
                %
                idxprev=1;
                for j=1:length(nodes),
                    node_elem{nodes(j)}=[node_elem{nodes(j)} idx(idxprev:idx2(j))];
                    idxprev=idx2(j)+1;
                end
            end
        end
        function [ids]=getBlockIDs(obj)
            ids=[obj.ID];
        end
        function [idx]=id2idx(obj,id)
            ids=[obj.ID];
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
        function nodes=NodesInBlocks(obj,blkidx)
            if nargin==1 || length(obj)==1,
                blkidx=1;
            end
            nnodes=zeros(length(blkidx),1);
            for i=1:length(blkidx),
                nnodes(i)=length(unique(obj(blkidx(i)).Connectivity(:)));
            end
            nodes=zeros(sum(nnodes),1);
            ii=1;
            for i=1:length(blkidx),
                nodes(ii:ii+nnodes(i)-1)=unique(obj(blkidx(i)).Connectivity(:));
                ii=ii+nnodes(i);
            end
            nodes=sort(nodes);
        end
        function var=get.NumElementsBlk(obj)
            var=zeros(length(obj),1);
            for i=1:length(obj),
                var(i)=size(obj.Connectivity,1);
            end
        end
    end
end
