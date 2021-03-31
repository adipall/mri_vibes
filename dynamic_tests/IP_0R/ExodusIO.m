classdef ExodusIO < handle
    %% A Exodus reader using matlab netcdf stuff
    
    % Author Todd Simmermacher
    % Tel: 844-0096
    % Last Updated: 12/15/2010
    % Sandia National Laboratories
    %
    properties
        Nid
        NumDimensions
        NumVariables
        NumGlobalAttributes
        UnlimitedDim
        ModelDim % 1d 2d or 3d
        DefaultMode ='nc_noclobber'  % change to nc_clobber if you like to overwrite files
        IO_WordSize  = 4 % 4, 8
        FileSize = 1  % 0 (old style) 1 large file size
        VarIDGlobal
        VarIDNodal
        VarIDElem
        NamesElem
        ElemTruthTable
        NodesetTruthTable
        SidesetTruthTable % Mike Added
        VarIDNS
        VarIDSS
        LastTime
    end
    properties (Access=private)
        Version       = 3.01
        APIVersion   = 4.01
        StrMaxLength = 33
        LineMaxLength = 81
        NemFileVersion = 2.6
        NemAPIVersion = 3.06
    end
    methods
        function obj=ExodusIO(filename,type)
            useNetcdf4 = true; 
            if nargin<2,
                type='rw';
            end
            switch type
                case 'rw',
                    obj.Nid=netcdf.open(filename,'write');
                    netcdf.setFill(obj.Nid,'nofill');  %%% I think this is correct.
                    obj=obj.updateProperties;
                    obj.ModelDim=obj.isDim('num_dim');
                    %% check to see how many dimensions the model has
                    noerror=1;
                    numdims=1;
                    while noerror,
                        try
                            netcdf.inqDim(obj.Nid,numdims);
                            numdims=numdims+1;
                        catch
                            noerror=0;
                        end
                    end
                    if numdims>=netcdf.getConstant('NC_MAX_DIMS'),
                        warning('Model may have too many features to write completely\n from Matlab: numdims = %d',numdims);
                    end
                case 'create'
                    %% default for ExodusIO is for 64bit_offset and noclobber
                    b64=netcdf.getConstant('64bit_offset');
                    nc=netcdf.getConstant(obj.DefaultMode);
                    %share=netcdf.getConstant('nc_share');
                    %%
% 					%%old style netcdf  --------------------------------
                    if ~useNetcdf4
                        obj.Nid=netcdf.create(filename,bitor(b64,nc));
                        % 					%-------------------------------------------------
                        
                    else
                        %netcdf4 option -------------- potential fix for 1024 dimensions issue found by Kevin Manktelow...needs some work still for compatibility
                        nccreate(filename,'netcdf4');
                        obj.Nid=netcdf.open(filename,'WRITE');
                        netcdf.reDef(obj.Nid);
                    end
					%-----------------------------------
                    netcdf.setFill(obj.Nid,'nofill');
                    %%mode = bitor();% bitwise or function
                    
                    %%%%
                    % add global attributes
                    %%%%
                    gbl=netcdf.getConstant('Global');
                    netcdf.putAtt(obj.Nid,gbl,'api_version',single(obj.APIVersion));
                    netcdf.putAtt(obj.Nid,gbl,'version',single(obj.Version));
                    if 1
                        netcdf.putAtt(obj.Nid,gbl,'nemesis_file_version',single(obj.NemFileVersion));
                        netcdf.putAtt(obj.Nid,gbl,'nemesis_api_version',single(obj.NemAPIVersion));
                    end
                    %%%
                    % Now add Dimensions
                    %%%
                    dimid=netcdf.defDim(obj.Nid,'len_string',obj.StrMaxLength);
                    dimid=netcdf.defDim(obj.Nid,'len_line',obj.LineMaxLength);
                    dimid=netcdf.defDim(obj.Nid,'four',4);
                   
                    netcdf.endDef(obj.Nid);
                otherwise
                    error('Unknown file open option')
            end
            obj=obj.updateProperties;
        end
        function delete(obj)
            %netcdf.close(obj.Nid);
        end
        function close(obj)
            netcdf.close(obj.Nid);
        end
        function [floating_wrd_sz,file_sz,ver,api_ver,lasttime]=rd_FileAttr(obj)
            %%
            %     [floating_wrd_sz,file_sz,ver,api_ver]=rd_FileAttr(obj)
            %%
            floating_wrd_sz=netcdf.getAtt(obj.Nid,netcdf.getConstant('NC_GLOBAL'),'floating_point_word_size');
            if obj.isAtt('file_size'),
                file_sz=netcdf.getAtt(obj.Nid,netcdf.getConstant('NC_GLOBAL'),'file_size');
            else
                file_sz=0;
            end
            obj.FileSize=file_sz;
            if obj.isAtt('last_written_time'),
                lasttime=netcdf.getAtt(obj.Nid,netcdf.getConstant('NC_GLOBAL'),'last_written_time');
            else
                lasttime=[];
            end
            obj.LastTime=lasttime;
            ver=netcdf.getAtt(obj.Nid,netcdf.getConstant('NC_GLOBAL'),'version');
            api_ver=netcdf.getAtt(obj.Nid,netcdf.getConstant('NC_GLOBAL'),'api_version');    
        end
        function title=rd_Title(obj)
            %%
            %   title=rd_Title(obj)
            %%
            title=netcdf.getAtt(obj.Nid,netcdf.getConstant('NC_GLOBAL'),'title');
        end
        function qa=rd_QA(obj)
            % to see QA records, squeeze(qa(:,:,i));
            [flag]=obj.isDim('num_qa_rec');
            if flag,
                vid=netcdf.inqVarID(obj.Nid,'qa_records');
                qa_hold=netcdf.getVar(obj.Nid,vid);
                [nls,n4,nqa]=size(qa_hold);
                qa=cell(1,nqa);
                for i=1:nqa,
                    qa{i}={deblank(squeeze(qa_hold(:,1,i))'),deblank(squeeze(qa_hold(:,2,i))'), ...
                        deblank(squeeze(qa_hold(:,3,i))'),deblank(squeeze(qa_hold(:,4,i))')};
                end
            else
                qa=[];
            end
        end
        function infos=rd_Info(obj)
            [flag,vid]=obj.isVar('info_records');
            infos=[];
            if flag,
                info=netcdf.getVar(obj.Nid,vid)';
                infos=cell(size(info,1),1);
                for i=1:size(info,1),
                    infos{i}=deblank(info(i,:));
                end
            end
        end
        function [coord,nmap,names]=rd_Nodes(obj)
            if obj.isAtt('file_size'),
                file_sz=netcdf.getAtt(obj.Nid,netcdf.getConstant('NC_GLOBAL'),'file_size');
            else
                file_sz=0;
            end
            if file_sz,
                idx=netcdf.inqVarID(obj.Nid,'coordx');
                if obj.ModelDim>1,
                    idy=netcdf.inqVarID(obj.Nid,'coordy');
                end
                if obj.ModelDim>2,
                    idz=netcdf.inqVarID(obj.Nid,'coordz');
                end
                switch obj.ModelDim
                    case 1
                        coord=netcdf.getVar(obj.Nid,idx);
                    case 2
                        coord=[netcdf.getVar(obj.Nid,idx) netcdf.getVar(obj.Nid,idy)];
                    case 3
                        coord=[netcdf.getVar(obj.Nid,idx) netcdf.getVar(obj.Nid,idy) netcdf.getVar(obj.Nid,idz)];
                    otherwise
                        error('Unknown Dimension: %d',obj.ModelDim)
                end
            else
                id=netcdf.inqVarID(obj.Nid,'coord');
                coord=netcdf.getVar(obj.Nid,id);
            end
            %%%
            [flag,vid]=obj.isVar('coor_names');
            name=netcdf.getVar(obj.Nid,vid);
            names=cell(obj.ModelDim,1);
            for i=1:obj.ModelDim,
                names{i}=deblank(name(:,i)');
            end
            %%%
            [flag,vid]=obj.isVar('node_num_map');
            nmap=[];
            if flag,
                nmap=netcdf.getVar(obj.Nid,vid);
            end
            
            %%% catch funny node ids from nemesis file (use original if
            %%% exist)
            [flag,vid]=obj.isVar('nmap_names');
            if flag
                %% see if name is original_global_id_map
                mapNameVal = netcdf.getVar(obj.Nid,vid);
                if length(mapNameVal) >= 22 && strcmp(mapNameVal(1:22)','original_global_id_map')
                    %% if that is it, try to get that node map and use it
                    [flag,vid] = obj.isVar('node_map1');
                    if flag
                        themap = netcdf.getVar(obj.Nid,vid);
                        
                        if length(themap) == length(nmap)
                            nmap = themap;
                        end
                    end
                end
            end
            
            
            
        end
        function blkids=rd_BlockIDs(obj)
            vid=netcdf.inqVarID(obj.Nid,'eb_prop1');
            blkids=netcdf.getVar(obj.Nid,vid);
        end
        function [conn,elemtype,name,status,att,attname]=rd_Block(obj,blkid)
            blkids=obj.rd_BlockIDs;
            idx=find(blkids==blkid);
            if isempty(idx),
                error('Block number %d not found',blkid);
            end
            %%%
            vid=netcdf.inqVarID(obj.Nid,'eb_status');
            stats=netcdf.getVar(obj.Nid,vid);
            status=stats(idx);
            %%
            [flag,vid]=obj.isVar(sprintf('connect%d',idx));
            if flag,
                conn=netcdf.getVar(obj.Nid,vid)';
                elemtype=netcdf.getAtt(obj.Nid,vid,'elem_type');
            else
                conn=[];
                elemtype=[];
            end
            %%
            [flag,vid]=obj.isVar('eb_names');
            if flag,
                names=netcdf.getVar(obj.Nid,vid)';
                name=names(idx,:)';
            else
                name=[];
            end
            [flag,vid]=obj.isVar(sprintf('attrib%d',idx));
            if flag,
                att=netcdf.getVar(obj.Nid,vid)';
            else
                att=[];
            end
            [flag,vid]=obj.isVar(sprintf('attrib_name%d',idx));
            if flag,
                [vname,xtype,dimids,nattr]=netcdf.inqVar(obj.Nid,vid);
                [dimname,dimlen]=netcdf.inqDim(obj.Nid,dimids(2));
                an=netcdf.getVar(obj.Nid,vid)';
                
                for i=1:dimlen,
                    attname{i}=an(i,:);
                end
            else
                attname=[];
            end
        end
        function emap=rd_EMap(obj)
            [flag,vid]=obj.isVar('elem_num_map');
            if flag,
                emap=netcdf.getVar(obj.Nid,vid);
            else
                emap=[];
            end
        end
        function nsids=rd_NodesetIDs(obj)
            [flag,vid]=obj.isVar('ns_prop1');
            if flag,
                nsids=netcdf.getVar(obj.Nid,vid);
            else
                nsids=[];
            end
        end
        function [nodes,dist_fact,name,status]=rd_Nodeset(obj,nsid)
            nsids=obj.rd_NodesetIDs;
            idx=find(nsids==nsid);
            if isempty(idx),
                error('Nodeset number %d not found',nsid);
            end
            %%%%%
            vid=netcdf.inqVarID(obj.Nid,'ns_status');
            stats=netcdf.getVar(obj.Nid,vid);
            status=stats(idx);
            %%%%%
            [flag,vid]=obj.isVar(sprintf('node_ns%d',idx));
            %vid=netcdf.inqVarID(obj.Nid,sprintf('node_ns%d',idx));
            if flag,
                nodes=netcdf.getVar(obj.Nid,vid);
            else
                nodes=[];
            end
            %%%%%
            [flag,vid]=obj.isVar(sprintf('dist_fact_ns%d',idx));
            if flag,
                dist_fact=netcdf.getVar(obj.Nid,vid);
            else
                dist_fact=[];
            end
            %%%%%
            [flag,vid]=obj.isVar('ns_names');
            if flag,
                names=netcdf.getVar(obj.Nid,vid)';
                name=names(idx,:);
            else
                name=[];
            end
        end
        function ssids=rd_SidesetIDs(obj)
            [flag,vid]=obj.isVar('ss_prop1');
            if flag,
                ssids=netcdf.getVar(obj.Nid,vid);
            else
                ssids=[];
            end
        end
        function [elems,sides,dist_fact,status,name]=rd_Sideset(obj,ssid)
            ssids=obj.rd_SidesetIDs;
            idx=find(ssid==ssids);
            if isempty(idx),
                error('Sideset number %d not found',ssid);
            end
            %%%%%
            vid=netcdf.inqVarID(obj.Nid,'ss_status');
            stats=netcdf.getVar(obj.Nid,vid);
            status=stats(idx);
            %%%%%
            [flag,vid]=obj.isVar(sprintf('elem_ss%d',idx));
            if flag,
                elems=netcdf.getVar(obj.Nid,vid);
            else
                elems=[];
            end
            %%%%%
            [flag,vid]=obj.isVar(sprintf('side_ss%d',idx));
            if flag,
                sides=netcdf.getVar(obj.Nid,vid);
            else
                sides=[];
            end
            %%%%%
            [flag,vid]=obj.isVar(sprintf('dist_fact_ss%d',idx));
            if flag,
                dist_fact=netcdf.getVar(obj.Nid,vid);
            else
                dist_fact=[];
            end
            %%%%%
            [flag,vid]=obj.isVar('ss_names');
            if flag,
                names=netcdf.getVar(obj.Nid,vid)';
                name=names(idx,:);
            else
                name=[];
            end
        end
        function [frame_struct]=rd_CoordFrames(obj)
            %%#define NUM_CFRAMES  "num_cframes"
            %%#define NUM_CFRAME9  "num_cframes_9"
            %%#define FRAME_COORDS "frame_coordinates"
            %%#define FRAME_IDS    "frame_ids"
            %%#define FRAME_TAGS   "frame_tags"
            [flag,vid]=obj.isVar('frame_ids');
            if flag,
                frame_ids=netcdf.getVar(obj.Nid,vid);
                varid=netcdf.inqVarID(obj.Nid,'frame_tags');
                frame_tags=netcdf.getVar(obj.Nid,varid);
                varid=netcdf.inqVarID(obj.Nid,'frame_coordinates');
                frame_coords=netcdf.getVar(obj.Nid,varid);
                ncoord=length(frame_ids);
                
                for i=1:ncoord,
                    frame_struct(i).id=frame_ids(i);
                    frame_struct(i).tag=deblank(frame_tags(i,:)');
                    frame_struct(i).coord_sys=reshape(frame_coords(9*i-8:9*i),3,3)';
                end
            else
                frame_struct=[];
            end
        end
        function nem_struct=rd_Nemesis(obj)
            datanames=obj.nem_datanames;
            dimnames=obj.nem_dimnames;
            nem_struct=[];
            
            getDebugInfo = 0;
            
            if getDebugInfo
                varIds = netcdf.inqVarIDs(obj.Nid);
                
                for ii = varIds
                    disp(netcdf.inqVar(obj.Nid,ii));
                    allData(ii+1).Name = netcdf.inqVar(obj.Nid,ii);
                    allData(ii+1).Data = netcdf.getVar(obj.Nid,ii);
                end
                
            end
            
            [flag]=isVar(obj,datanames{1});
            if flag
                for i=1:length(datanames),
                    if ~isempty(datanames{i}),
                        try
                            vid=netcdf.inqVarID(obj.Nid,datanames{i});
                            nem_struct.(datanames{i})=netcdf.getVar(obj.Nid,vid);
                        end
                    else
                        try
                            nem_struct.(dimnames{i})=obj.isDim(dimnames{i});
                        end
                    end
                end
            end
        end
        function time=rd_Time(obj)
            vid=netcdf.inqVarID(obj.Nid,'time_whole');
            time=netcdf.getVar(obj.Nid,vid);
        end
        function gnames=rd_GlobalVarsNames(obj)
            [flag,vid]=obj.isVar('name_glo_var');
            if flag,
                gn=netcdf.getVar(obj.Nid,vid)';
                gnames=cell(1,size(gn,1));
                for i=1:size(gn,1),
                  gnames{i}=deblank(gn(i,:));
                end
            else
                gnames=[];
            end
        end
        function data=rd_GlobalVar(obj,name)
            gnames=obj.rd_GlobalVarsNames;
            idx=strmatch(name,gnames,'exact');
            if isempty(idx),
                data=[];
                return
            end
            vid=netcdf.inqVarID(obj.Nid,'vals_glo_var');
            hold=netcdf.getVar(obj.Nid,vid);
            data=squeeze(hold(idx,:));
        end
        function data=rd_GlobalVarByInd(obj,idx)
%             gnames=obj.rd_GlobalVarsNames;
%             idx=strmatch(name,gnames,'exact');
            if isempty(idx),
                data=[];
                return
            end
            vid=netcdf.inqVarID(obj.Nid,'vals_glo_var');
            hold=netcdf.getVar(obj.Nid,vid);
            data=squeeze(hold(idx,:));
        end
        function nnames=rd_NodalVarsNames(obj)
            [flag,vid]=obj.isVar('name_nod_var');
            if flag,
                nn=netcdf.getVar(obj.Nid,vid)';
                nnames=cell(1,size(nn,1));
                for i=1:size(nn,1)
                    nnames{i}=deblank(nn(i,:));
                end
            else
                nnames=[];
            end
        end
        function [data]=rd_NodalVar(obj,name)
            nnames=obj.rd_NodalVarsNames;
            idx=strmatch(name,nnames,'exact');
            if isempty(idx),
                data=[];
                return
            end
            
            if obj.FileSize,
                vid=netcdf.inqVarID(obj.Nid,sprintf('vals_nod_var%d',idx));
                data=netcdf.getVar(obj.Nid,vid);
            else
                vid=netcdf.inqVarID(obj.Nid,'vals_nod_var');
                hold=netcdf.getVar(obj.Nid,vid);
                data=squeeze(hold(:,idx,:));
            end
            
        end
        function enames=rd_ElemVarsNames(obj)
            [flag,vid]=obj.isVar('name_elem_var');
            if flag,
                en=netcdf.getVar(obj.Nid,vid)';
                enames=cell(1,size(en,1));
                for i=1:size(en,1),
                    enames{i}=deblank(en(i,:));
                end
            else
                enames=[];
            end
        end
        function [data]=rd_ElemVar(obj,name,blkid)
            enames=obj.rd_ElemVarsNames;
            idxn=strmatch(name,enames,'exact');
            if isempty(idxn),
                error('Element Variable %s not found',name);
            end
            blkids=obj.rd_BlockIDs;
            idxb=find(blkids==blkid);
            if isempty(idxb),
                error('Block number %d not found',blkid);
            end
            [flag,vid]=obj.isVar(sprintf('vals_elem_var%deb%d',idxn,idxb));
            if flag
                data=netcdf.getVar(obj.Nid,vid);
            else
                data=[];
            end
        end
        function nsnames=rd_NodesetVarsNames(obj)
            [flag,vid]=obj.isVar('name_nset_var');
            if flag,
                nsn=netcdf.getVar(obj.Nid,vid)';
                nsnames=cell(1,size(nsn,1));
                for i=1:size(nsn,1),
                    nsnames{i}=deblank(nsn(i,:));
                end
            else
                nsnames=[];
            end
        end
        function [data]=rd_NodesetVar(obj,name,nsid)
            nsnames=obj.rd_NodesetVarsNames;
            idxn=strmatch(name,nsnames,'exact');
            if isempty(idxn),
                error('Nodeset Variable %s not found',name);
            end
            nsids=obj.rd_NodesetIDs;
            idxns=find(nsids==nsid);
            if isempty(idxns),
                error('Nodeset number %d not found',nsid);
            end
            [flag,vid]=obj.isVar(sprintf('vals_nset_var%dns%d',idxn,idxns));
            if flag
                data=netcdf.getVar(obj.Nid,vid);
            else
                data=[];
            end
        end
        %------------------------------------------------------------------
        % Mike Added
        function ssnames=rd_SidesetVarsNames(obj)
            [flag,vid]=obj.isVar('name_sset_var');
            if flag,
                ssn=netcdf.getVar(obj.Nid,vid)';
                ssnames=cell(1,size(ssn,1));
                for i=1:size(ssn,1),
                    ssnames{i}=deblank(ssn(i,:));
                end
            else
                ssnames=[];
            end
        end
        function [data]=rd_SidesetVar(obj,name,ssid)
            ssnames=obj.rd_SidesetVarsNames;
            idxs=strmatch(name,ssnames,'exact');
            if isempty(idxs),
                error('Sideset Variable %s not found',name);
            end
            ssids=obj.rd_SidesetIDs;
            idxss=find(ssids==ssid);
            if isempty(idxss),
                error('Sideset number %d not found',ssid);
            end
            [flag,vid]=obj.isVar(sprintf('vals_sset_var%dss%d',idxs,idxss));
            if flag
                data=netcdf.getVar(obj.Nid,vid);
            else
                data=[];
            end
        end
        
        %------------------------------------------------------------------
        %%%%
        function [Kr,Cr,Mr,cbmap,otm,outmap,otme,outemap, ...
                ddof,din,deig,dnc]=rd_SuperElement(obj)
            [flag,vid]=obj.isVar('Kr');
            if flag,
                Kr=netcdf.getVar(obj.Nid,vid)';
            else
                Kr=[];
            end
            [flag,vid]=obj.isVar('Cr');
            if flag,
                Cr=netcdf.getVar(obj.Nid,vid)';
            else
                Cr=[];
            end
            [flag,vid]=obj.isVar('Mr');
            if flag,
                Mr=netcdf.getVar(obj.Nid,vid)';
            else
                Mr=[];
            end
            [flag,vid]=obj.isVar('cbmap');
            if flag,
                cbmap=netcdf.getVar(obj.Nid,vid)';
            else
                cbmap=[];
            end
            [flag,vid]=obj.isVar('OTM');
            if flag,
                otm=netcdf.getVar(obj.Nid,vid)';
            else
                otm=[];
            end
            [flag,vid]=obj.isVar('OutMap');
            if flag,
                outmap=netcdf.getVar(obj.Nid,vid);
            else
                outmap=[];
            end
            [flag,vid]=obj.isVar('OTME');
            if flag,
                otme=netcdf.getVar(obj.Nid,vid)';
            else
                otme=[];
            end
            [flag,vid]=obj.isVar('OutElemMap');
            if flag,
                outemap=netcdf.getVar(obj.Nid,vid);
            else
                outemap=[];
            end
            ddof=obj.isDim('NumDof');
            din=obj.isDim('NumInterfaceNodes');
            deig=obj.isDim('NumEig');
            if isempty(deig),
                deig=0;
            end
            dnc=obj.isDim('NumConstraints');
        end
        %%%%%
        function obj=create(obj,io_ws,file_size,lasttime, ...
                title,num_dim,num_nodes,num_elem,num_elem_blk, ...
                num_nodesets,num_sidesets)
            netcdf.reDef(obj.Nid);
            %%
            gbl=netcdf.getConstant('Global');
            
            netcdf.putAtt(obj.Nid,gbl,'floating_point_word_size',int32(io_ws));
            obj.IO_WordSize=io_ws;
            netcdf.putAtt(obj.Nid,gbl,'file_size',int32(file_size));
            obj.FileSize=file_size;
            dtime=netcdf.defDim(obj.Nid,'time_step',netcdf.getConstant('NC_UNLIMITED'));
            [ju,type]=obj.setType(1);
            netcdf.defVar(obj.Nid,'time_whole',type,dtime);
            netcdf.putAtt(obj.Nid,netcdf.getConstant('Global'),'title',title');
            netcdf.putAtt(obj.Nid,netcdf.getConstant('Global'),'number_equations',int32(0));
            if ~isempty(lasttime),
                netcdf.putAtt(obj.Nid,netcdf.getConstant('Global'),'last_written_time',lasttime);
            end
            dimid=netcdf.defDim(obj.Nid,'num_dim',num_dim);
            dimid=netcdf.defDim(obj.Nid,'num_nodes',num_nodes);
            dimid=netcdf.defDim(obj.Nid,'num_elem',num_elem);
            dimid=netcdf.defDim(obj.Nid,'num_el_blk',num_elem_blk);
            if num_nodesets,
                dimid=netcdf.defDim(obj.Nid,'num_node_sets',num_nodesets);
            end
            if num_sidesets,
                dimid=netcdf.defDim(obj.Nid,'num_side_sets',num_sidesets);
            end
            % end define mode
            netcdf.endDef(obj.Nid);
            % update number of dimensions
            obj=obj.updateProperties;
        end
        function obj=wr_QA(obj,qa_records)
            %
            % qa{1}={'qaline1','qaline2','qaline3','qaline4'};
            % qa{2}={'test1','test2','test3','test4'};
            %
            
            % reorganize qa
            nqa=length(qa_records);
            qa=repmat(char(0),[obj.StrMaxLength,4,nqa]);
            for i=1:nqa,
                if ~isempty(qa_records{i}{1}),
                    qa(1:length(qa_records{i}{1}),1,i)=qa_records{i}{1}';
                end
                if ~isempty(qa_records{i}{2}),
                    qa(1:length(qa_records{i}{2}),2,i)=qa_records{i}{2}';
                end
                if ~isempty(qa_records{i}{3}),
                    qa(1:length(qa_records{i}{3}),3,i)=qa_records{i}{3}';
                end
                if ~isempty(qa_records{i}{4}),
                    qa(1:length(qa_records{i}{4}),4,i)=qa_records{i}{4}';
                end
            end
            % Go into define mode
            netcdf.reDef(obj.Nid);
            dim=isDim(obj,'num_qa_rec');
            if ~isempty(dim),
                error('QA records already exist');
            end
            nqadim=netcdf.defDim(obj.Nid,'num_qa_rec',nqa);
            %
            [a,n4dim]=isDim(obj,'four');
            [a,nlsdim]=isDim(obj,'len_string');
            %%
            vid=netcdf.defVar(obj.Nid,'qa_records','char',[nlsdim,n4dim,nqadim]);
            %
            % exit define mode
            %
            netcdf.endDef(obj.Nid);
            %%
            netcdf.putVar(obj.Nid,vid,qa);
            obj=obj.updateProperties;
        end
        function obj=wr_Info(obj,info_records)
            %
            % info={'infoline1','infoline2','infoline3'};
            %
            % reorganize info
            ninfo=length(info_records);
%             info=repmat(' ',[obj.LineMaxLength,ninfo]);
%             for i=1:ninfo,
%                 info(1:length(info_records{i}),i)=info_records{i};
%             end
            info=ExodusIO.mk_str_array(info_records,obj.LineMaxLength);
            % Go into define mode
            netcdf.reDef(obj.Nid);
            dim=isDim(obj,'num_info');
            if ~isempty(dim),
                error('Info records already exist');
            end
            ninfodim=netcdf.defDim(obj.Nid,'num_info',ninfo);
            %
            [a,nlldim]=isDim(obj,'len_line');
            %%
            vid=netcdf.defVar(obj.Nid,'info_records','char',[nlldim,ninfodim]);
            %
            % exit define mode
            %
            netcdf.endDef(obj.Nid);
            %%
            netcdf.putVar(obj.Nid,vid,info);
            obj=obj.updateProperties;
        end
        function obj=wr_Nodes(obj,coord,nmap,names)
            %
            %
            %
            [ls,nlsdim]=isDim(obj,'len_string');
            [dnodes,dnodid]=obj.isDim('num_nodes');
            [ddim,ddimid]=obj.isDim('num_dim');

            netcdf.reDef(obj.Nid);
            
            vinames=netcdf.defVar(obj.Nid,'coor_names','char',[nlsdim,ddimid]);
            if nargin<2 || isempty(nmap),
            else
                vidmap=netcdf.defVar(obj.Nid,'node_num_map','int',dnodid);
            end
            [coord,type]=obj.setType(coord);
            if obj.FileSize,
                vidx=netcdf.defVar(obj.Nid,'coordx',type,dnodid);
                if ddim>1,
                    vidy=netcdf.defVar(obj.Nid,'coordy',type,dnodid);
                end
                if ddim>2,
                    vidz=netcdf.defVar(obj.Nid,'coordz',type,dnodid);
                end
            else
                vidx=netcdf.defVar(obj.Nid,'coord',type,[dnodid,ddimid]);
            end
            % end define mode
            netcdf.endDef(obj.Nid);
            %
           
            if obj.FileSize,
                netcdf.putVar(obj.Nid,vidx,coord(:,1));
                if ddim>1,
                    netcdf.putVar(obj.Nid,vidy,coord(:,2));
                end
                if ddim>2,
                    netcdf.putVar(obj.Nid,vidz,coord(:,3));
                end
            else
                netcdf.putVar(obj.Nid,vidx,coord);
            end
            if nargin<2 || isempty(nmap),
            else
                netcdf.putVar(obj.Nid,vidmap,int32(nmap));
            end
            att=ExodusIO.mk_str_array(names,obj.StrMaxLength);
            netcdf.putVar(obj.Nid,vinames,att)
            %% update file info
            obj=obj.updateProperties;
        end
        function obj=wr_Blocks(obj,blkids,conn,elemtype,blkname,status,att,attname)
            % blkids = nblks x 1 matrix
            % conn = nblks cell array , each cell nelem per blk x nnodes per block
            % elemtyp = cell array of strings
            % name = cell array of strings or empty
            % status = nblks x 1 matrix (1 or zero)
            % att = cell array of attributes (each cell has a matrix for
            % each block
            % attname = cell array name of attribute
            
            netcdf.reDef(obj.Nid);
            [a,nlsdim]=isDim(obj,'len_string');
            %% define dimensions and variables
            [dblks,diblks]=obj.isDim('num_el_blk');
            if isempty(dblks),
                error('File was improperly initialized')
            end
            if ~isempty(status),
                vistat=netcdf.defVar(obj.Nid,'eb_status','int',diblks);
            end
            if ~isempty(blkname),
                vblkn=netcdf.defVar(obj.Nid,'eb_names','char',[nlsdim,diblks]);
            end
            viblks=netcdf.defVar(obj.Nid,'eb_prop1','int',diblks);
            netcdf.putAtt(obj.Nid,viblks,'name','ID');
            diel=zeros(length(blkids),1);
            dinod=zeros(length(blkids),1);
            diatt=zeros(length(blkids),1);
            viblk=zeros(length(blkids),1);
            viatt=zeros(length(blkids),1);
            viattn=zeros(length(blkids),1);
            %
            [jnk,type]=obj.setType(1);
            %
            for i=1:length(blkids)
                % define dimensions and variables
                if ~isempty(conn{i}),
                    diel(i)=netcdf.defDim(obj.Nid,sprintf('num_el_in_blk%d',i),size(conn{i},1));
                    dinod(i)=netcdf.defDim(obj.Nid,sprintf('num_nod_per_el%d',i),size(conn{i},2));
                    viblk(i)=netcdf.defVar(obj.Nid,sprintf('connect%d',i),'int',[dinod(i),diel(i)]);
                    netcdf.putAtt(obj.Nid,viblk(i),'elem_type',elemtype{i});
                    %%
                    if nargin>6 && ~isempty(att) && ~isempty(att{i}),
                        diatt(i)=netcdf.defDim(obj.Nid,sprintf('num_att_in_blk%d',i),size(att{i},2));
                        viatt(i)=netcdf.defVar(obj.Nid,sprintf('attrib%d',i),type,[diatt(i),diel(i)]);
                        if nargin>7 && ~isempty(attname) && ~isempty(attname{i}),
                            viattn(i)=netcdf.defVar(obj.Nid,sprintf('attrib_name%d',i),'char',[nlsdim,diatt(i)]);
                        end
                    end
                end
            end
            %% Populate variables
            netcdf.endDef(obj.Nid);
            %%%
            if ~isempty(status),
                netcdf.putVar(obj.Nid,vistat,status);
            end
            netcdf.putVar(obj.Nid,viblks,int32(blkids));
            %
            if ~isempty(blkname),
                arr_blkname=ExodusIO.mk_str_array(blkname,obj.StrMaxLength);
                netcdf.putVar(obj.Nid,vblkn,arr_blkname);
            end
            %
            
            for i=1:length(blkids),
                if ~isempty(conn{i}),
                    netcdf.putVar(obj.Nid,viblk(i),int32(conn{i}'));
                    
                    if nargin>6 && ~isempty(att) && ~isempty(att{i}),
                        
                        netcdf.putVar(obj.Nid,viatt(i),att{i}');
                        if nargin>7 && ~isempty(attname) && ~isempty(attname{i})
                            atn=ExodusIO.mk_str_array(attname{i},obj.StrMaxLength);
                            netcdf.putVar(obj.Nid,viattn(i),atn);
                        end
                    end
                end
                
            end
            obj=obj.updateProperties;
        end
        function obj=wr_EMap(obj,emap)
            [dlen,dim]=obj.isDim('num_elem');
            netcdf.reDef(obj.Nid);
            vid=netcdf.defVar(obj.Nid,'elem_num_map','int',dim);
            netcdf.endDef(obj.Nid);
            netcdf.putVar(obj.Nid,vid,emap);
        end
        function obj=wr_Nodesets(obj,ids,nodes,dist_fact,name,status)
            % ids = nsnum x 1 matrix
            % nodes = nns cell array , each cell nnode per ns long
            % dist_fact = cell array of strings
            % name = cell array of strings or empty
            % status = numns x 1 matrix (1 or zero)
            
            
            netcdf.reDef(obj.Nid);
            [a,nlsdim]=isDim(obj,'len_string');
            %% define dimensions and variables
            [dns,dins]=obj.isDim('num_node_sets');
            if isempty(dns),
                error('File was improperly initialized')
            end
            if ~isempty(status),
                vistat=netcdf.defVar(obj.Nid,'ns_status','int',dins);
            end
            if ~isempty(name),
                vnsn=netcdf.defVar(obj.Nid,'ns_names','char',[nlsdim,dins]);
            end
            vinsp=netcdf.defVar(obj.Nid,'ns_prop1','int',dins);
            netcdf.putAtt(obj.Nid,vinsp,'name','ID');
            dinod=zeros(length(ids),1);
            vins=zeros(length(ids),1);
            vidf=zeros(length(ids),1);
            %
            [jnk,type]=obj.setType(1);
            %
            for i=1:length(ids)
                % define dimensions and variables
                if ~isempty(nodes{i}),
                    dinod(i)=netcdf.defDim(obj.Nid,sprintf('num_nod_ns%d',i),length(nodes{i}));
                    vins(i)=netcdf.defVar(obj.Nid,sprintf('node_ns%d',i),'int',dinod(i));
                    vidf(i)=netcdf.defVar(obj.Nid,sprintf('dist_fact_ns%d',i),type,dinod(i));
                end
                %%
            end
            %% Populate variables
            netcdf.endDef(obj.Nid);
            %%%
            if ~isempty(status),
                netcdf.putVar(obj.Nid,vistat,status);
            end
            if ~isempty(name),
                arr_name=ExodusIO.mk_str_array(name,obj.StrMaxLength);
                netcdf.putVar(obj.Nid,vnsn,arr_name);
            end
            netcdf.putVar(obj.Nid,vinsp,int32(ids));
            %
            for i=1:length(ids),
                if ~isempty(nodes{i}),
                    netcdf.putVar(obj.Nid,vins(i),int32(nodes{i}));
                    if ~isempty(dist_fact{i}),
                        netcdf.putVar(obj.Nid,vidf(i),dist_fact{i});
                    end
                end
            end
            obj=obj.updateProperties;
        end
        function obj=wr_Sidesets(obj,ids,elems,sides,dist_fact,name,status)
            % ids = nsnum x 1 matrix
            % elems = nssets cell array , each cell nelemss per ss long
            % sides = nssets cell array , each cell nelemss per ss long
            % dist_fact = cell array of strings
            % name = cell array of strings or empty
            % status = numns x 1 matrix (1 or zero) 
            
            netcdf.reDef(obj.Nid);
            [a,nlsdim]=isDim(obj,'len_string');
            %% define dimensions and variables
            [dss,diss]=obj.isDim('num_side_sets');
            if isempty(diss),
                error('File was improperly initialized')
            end
            if ~isempty(status),
                vistat=netcdf.defVar(obj.Nid,'ss_status','int',diss);
            end
            if ~isempty(name),
                vssn=netcdf.defVar(obj.Nid,'ss_names','char',[nlsdim,diss]);
            end
            vissp=netcdf.defVar(obj.Nid,'ss_prop1','int',diss);
            netcdf.putAtt(obj.Nid,vissp,'name','ID');
            diel=zeros(length(ids),1);
            didf=zeros(length(ids),1);
            visse=zeros(length(ids),1);
            visss=zeros(length(ids),1);
            vidf=zeros(length(ids),1);
            %
            [jnk,type]=obj.setType(1);
            %
            for i=1:length(ids)
                % define dimensions and variables
                if ~isempty(elems{i}),
                    diel(i)=netcdf.defDim(obj.Nid,sprintf('num_side_ss%d',i),length(elems{i}));
                    visse(i)=netcdf.defVar(obj.Nid,sprintf('elem_ss%d',i),'int',diel(i));
                    visss(i)=netcdf.defVar(obj.Nid,sprintf('side_ss%d',i),'int',diel(i));
                    if nargin>4 && ~isempty(dist_fact) && ~isempty(dist_fact{i}),
                        didf(i)=netcdf.defDim(obj.Nid,sprintf('num_df_ss%d',i),length(dist_fact{i}));
                        vidf(i)=netcdf.defVar(obj.Nid,sprintf('dist_fact_ss%d',i),type,didf(i));
                    end
                end
                %%
            end
            %% Populate variables
            netcdf.endDef(obj.Nid);
            %%%
            if ~isempty(status),
                netcdf.putVar(obj.Nid,vistat,status);
            end
            if ~isempty(name),
                arr_name=ExodusIO.mk_str_array(name,obj.StrMaxLength);
                netcdf.putVar(obj.Nid,vssn,arr_name);
            end
            netcdf.putVar(obj.Nid,vissp,int32(ids));
            %
            for i=1:length(ids),
                if ~isempty(elems{i}),
                    netcdf.putVar(obj.Nid,visse(i),int32(elems{i}));
                    netcdf.putVar(obj.Nid,visss(i),int32(sides{i}));
                    if nargin>4 && ~isempty(dist_fact) && ~isempty(dist_fact{i}),
                        netcdf.putVar(obj.Nid,vidf(i),dist_fact{i});
                    end
                end
            end
            obj=obj.updateProperties;
        end
        function obj=wr_CoordFrames(obj,frame_struct)
             %%#define NUM_CFRAMES  "num_cframes"
            %%#define NUM_CFRAME9  "num_cframes_9"
            %%#define FRAME_COORDS "frame_coordinates"
            %%#define FRAME_IDS    "frame_ids"
            %%#define FRAME_TAGS   "frame_tags"
            netcdf.reDef(obj.Nid);
            [ju,type]=obj.setType(1);
            [ls,nlsdim]=isDim(obj,'len_string');
            dimf=netcdf.defDim(obj.Nid,'num_cframes',length(frame_struct));
            dim9=netcdf.defDim(obj.Nid,'num_cframes_9',9*length(frame_struct));
            varc=netcdf.defVar(obj.Nid,'frame_coordinates',type,dim9);
            vart=netcdf.defVar(obj.Nid,'frame_tags','char',[dimf,nlsdim]);
            varid=netcdf.defVar(obj.Nid,'frame_ids','int',dimf);
            
            netcdf.endDef(obj.Nid);
            frames=zeros(9*length(frame_struct),1);
            ids=zeros(length(frame_struct),1);
            cells=cell(length(frame_struct),1);
            for i=1:length(frame_struct),
                a=frames_struct(i).coord_sys';
                frames(9*i-8:9*i-1)=a(:);
                ids(i)=frame_struct(i).ids;
                cells{i}=frames_struct(i).tag;
            end
            netcdf.putVar(obj.Nid,varc,frames);
            netcdf.putVar(obj.Nid,varid,ids);
            arr=ExodusIO.mk_str_array(cells,obj.StrMaxLength);
            netcdf.putVar(obj.Nid,vart,arr);
            obj=obj.updateProperties;
        end
        function obj=wr_Nemesis(obj,nem_struct)
            nem_dim=obj.nem_dimnames;
            nem_data=obj.nem_datanames;
            
            netcdf.reDef(obj.Nid);
            var=zeros(length(nem_data),1);
            dim=zeros(length(nem_data),1);
            for i=1:length(nem_data),
                if isfield(nem_struct,nem_data{i}),
                    if isempty(nem_dim{i}),
                        var(i)=netcdf.defVar(obj.Nid,nem_data{i},'int',[]);
                        continue
                    end
                    idx=strmatch(nem_dim{i},nem_dim(1:i-1),'exact');
                    if i==1 || isempty(idx),
                        if ~isempty(nem_data{i}),
                            dim(i)=netcdf.defDim(obj.Nid,nem_dim{i},length(nem_struct.(nem_data{i})));
                        else
                            dim(i)=netcdf.defDim(obj.Nid,nem_dim{i},nem_struct.(nem_dim{i}));
                        end
                    else
                        dim(i)=dim(idx(1));
                    end
                    if ~isempty(nem_data{i}),
                        var(i)=netcdf.defVar(obj.Nid,nem_data{i},'int',dim(i));
                    end
                end
            end
            
            netcdf.endDef(obj.Nid);
            
            for i=1:length(nem_data),
                if isfield(nem_struct,nem_data{i}) && ~isempty(nem_data{i}) ,
                    netcdf.putVar(obj.Nid,var(i),nem_struct.(nem_data{i}));
                else
                    
                end
            end
            obj=obj.updateProperties;
        end
        function obj=wr_Time(obj,time)
            if isempty(time),
                return;
            end
            [flag,varid]=obj.isVar('time_whole');
            [dim]=obj.isDim('time_step');
            if flag,
                netcdf.putVar(obj.Nid,varid,dim,length(time),time);
            else
                error('Time step is not defined in Exodus File')
            end
            
            obj=obj.updateProperties;
        end
        function obj=init_GlobalVars(obj,names)
            [ls,nlsdim]=isDim(obj,'len_string');
            [ju,type]=obj.setType(1);
            
            netcdf.reDef(obj.Nid);
            dglo=netcdf.defDim(obj.Nid,'num_glo_var',length(names));
            vn=netcdf.defVar(obj.Nid,'name_glo_var','char',[nlsdim,dglo]);
            [ju,dt]=isDim(obj,'time_step');
            obj.VarIDGlobal=netcdf.defVar(obj.Nid,'vals_glo_var',type,[dglo,dt]);
            netcdf.endDef(obj.Nid);
            att=ExodusIO.mk_str_array(names,obj.StrMaxLength);
            netcdf.putVar(obj.Nid,vn,att);
        end
        function obj=init_NodalVars(obj,names)
            [ls,nlsdim]=isDim(obj,'len_string');
            [ju,type]=obj.setType(1);
            
            netcdf.reDef(obj.Nid);
            dno=netcdf.defDim(obj.Nid,'num_nod_var',length(names));

            vn=netcdf.defVar(obj.Nid,'name_nod_var','char',[nlsdim,dno]);
            [ju,dnnodes]=obj.isDim('num_nodes');
          
            [ju,dt]=isDim(obj,'time_step');
            
            if obj.FileSize,
                obj.VarIDNodal=zeros(length(names),1);
                for i=1:length(names),
                    obj.VarIDNodal(i)=netcdf.defVar(obj.Nid,sprintf('vals_nod_var%d',i),type,[dnnodes,dt]);
                end
            else
                obj.VarIDNodal=netcdf.defVar(obj.Nid,'vals_nod_var',type,[dnnodes,dno,dt]);
            end
            netcdf.endDef(obj.Nid);
            att=ExodusIO.mk_str_array(names,obj.StrMaxLength);
            netcdf.putVar(obj.Nid,vn,att);
        end
        function obj=init_ElemVars(obj,names,var_table)
            %% var_table is length(blkids)xlength(names) all zeros except
            %% in the locations where the blkid has the appropriate
            %% variable name
            [ls,nlsdim]=isDim(obj,'len_string');
            [ju,type]=obj.setType(1);
            obj.ElemTruthTable=var_table;
            
            netcdf.reDef(obj.Nid);
            [nblks,dblks]=obj.isDim('num_el_blk');
            del=netcdf.defDim(obj.Nid,'num_elem_var',length(names));
            vn=netcdf.defVar(obj.Nid,'name_elem_var','char',[nlsdim,del]);
            [ju,dt]=isDim(obj,'time_step');
            obj.VarIDElem=-ones(length(names),nblks);
            for i=1:nblks,
                for j=1:length(names),
                    if obj.ElemTruthTable(j,i),
                        [ju,delblk]=obj.isDim(sprintf('num_el_in_blk%d',i));
                        obj.VarIDElem(j,i)=netcdf.defVar(obj.Nid,sprintf('vals_elem_var%deb%d',j,i),type,[delblk,dt]);
                    end
                end
            end
            vvt=netcdf.defVar(obj.Nid,'elem_var_tab','int',[del,dblks]);
            netcdf.endDef(obj.Nid);
            att=ExodusIO.mk_str_array(names,obj.StrMaxLength);
            netcdf.putVar(obj.Nid,vn,att);
            netcdf.putVar(obj.Nid,vvt,var_table);
        end
        function obj=init_NodesetVars(obj,names,var_table)
            %% var_table is length(blkids)xlength(names) all zeros except
            %% in the locations where the blkid has the appropriate
            %% variable name
            [ls,nlsdim]=isDim(obj,'len_string');
            [ju,type]=obj.setType(1);
            obj.NodesetTruthTable=var_table;
            
            netcdf.reDef(obj.Nid);
            [nnsets,dns]=obj.isDim('num_node_sets');
            del=netcdf.defDim(obj.Nid,'num_nset_var',length(names));
            vn=netcdf.defVar(obj.Nid,'name_nset_var','char',[nlsdim,del]);
            [ju,dt]=isDim(obj,'time_step');
            obj.VarIDNS=-ones(length(names),nnsets);
            for i=1:nnsets,
                for j=1:length(names),
                    if obj.NodesetTruthTable(j,i),
                        [ju,delns]=obj.isDim(sprintf('num_nod_ns%d',i));
                        obj.VarIDNS(j,i)=netcdf.defVar(obj.Nid,sprintf('vals_nset_var%dns%d',j,i),type,[delns,dt]);
                    end
                end
            end
            vvt=netcdf.defVar(obj.Nid,'nset_var_tab','int',[del,dns]);
            netcdf.endDef(obj.Nid);
            att=ExodusIO.mk_str_array(names,obj.StrMaxLength);
            netcdf.putVar(obj.Nid,vn,att);
            netcdf.putVar(obj.Nid,vvt,var_table);
        end
        
        %------------------------------------------------------------------
        % Mike Added
        function obj=init_SidesetVars(obj,names,var_table)
            %% var_table is length(blkids)xlength(names) all zeros except
            %% in the locations where the blkid has the appropriate
            %% variable name
            [ls,nlsdim]=isDim(obj,'len_string');
            [ju,type]=obj.setType(1);
            obj.SidesetTruthTable=var_table;
            
            netcdf.reDef(obj.Nid);
            [nssets,dss]=obj.isDim('num_side_sets');
            del=netcdf.defDim(obj.Nid,'num_sset_var',length(names));
            vn=netcdf.defVar(obj.Nid,'name_sset_var','char',[nlsdim,del]);
            [ju,dt]=isDim(obj,'time_step');
            obj.VarIDSS=-ones(length(names),nssets);
            for i=1:nssets,
                for j=1:length(names),
                    if obj.SidesetTruthTable(j,i),
                        [ju,delns]=obj.isDim(sprintf('num_side_ss%d',i));
                        obj.VarIDSS(j,i)=netcdf.defVar(obj.Nid,sprintf('vals_sset_var%dss%d',j,i),type,[delns,dt]);
                    end
                end
            end
            vvt=netcdf.defVar(obj.Nid,'sset_var_tab','int',[del,dss]);
            netcdf.endDef(obj.Nid);
            att=ExodusIO.mk_str_array(names,obj.StrMaxLength);
            netcdf.putVar(obj.Nid,vn,att);
            netcdf.putVar(obj.Nid,vvt,var_table);
        end       
        %------------------------------------------------------------------
        
        function obj=wr_GlobalVars(obj,data)
            %netcdf.putVar(obj.Nid,obj.VarIDGlobal,time_steps(1),length(time_steps),data);
            
            netcdf.putVar(obj.Nid,obj.VarIDGlobal,data);
        end
        function obj=wr_NodalVars(obj,data)
            [ju1,nvars,ju2]=size(data);
            if obj.FileSize,
                for i=1:nvars,
                    %netcdf.putVar(obj.Nid,obj.VarIDNodal(i),1,length(time_steps),squeeze(data(:,i,:)));
                    netcdf.putVar(obj.Nid,obj.VarIDNodal(i),squeeze(data(:,i,:)));
                end
            else
                %netcdf.putVar(obj.Nid,obj.VarIDNodal,time_steps(1),length(time_steps),data);
                netcdf.putVar(obj.Nid,obj.VarIDNodal,data);
            end
        end
        function obj=wr_ElemVars(obj,data)
            [nvar,nblks]=size(data);
            for i=1:nblks,
                for j=1:nvar,
                    if obj.ElemTruthTable(j,i),
                        netcdf.putVar(obj.Nid,obj.VarIDElem(j,i),data{j,i});
                    end
                end
            end
        end
        function obj=wr_NodesetVars(obj,data)
            [nvar,nns]=size(data);
            for i=1:nns,
                for j=1:nvar,
                    if obj.NodesetTruthTable(j,i),
                        netcdf.putVar(obj.Nid,obj.VarIDNS(j,i),data{j,i});
                    end
                end
            end
        end
        %------------------------------------------------------------------
        % Mike Added
        function obj=wr_SidesetVars(obj,data)
            [svar,nss]=size(data);
            for i=1:nss,
                for j=1:svar,
                    if obj.SidesetTruthTable(j,i),
                        netcdf.putVar(obj.Nid,obj.VarIDSS(j,i),data{j,i});
                    end
                end
            end
        end
        %------------------------------------------------------------------
        
        %%%%
        function obj=wr_SuperElement(obj,Kr,Cr,Mr,cbmap,otm,outmap,otme, ...
                outemap,din,deig,dnc)
            netcdf.reDef(obj.Nid);
            if isempty(cbmap),
                return;
            end
            if ~isempty(Kr),
                n=length(Kr);
            elseif ~isempty(Mr),
                n=length(Mr);
            else
                error('Both Mr and Kr are empty')
            end
            dim_dof=netcdf.defDim(obj.Nid,'NumDof',n);
            dim_2=netcdf.defDim(obj.Nid,'two',2);
            if ~isempty(din),
                dim_node=netcdf.defDim(obj.Nid,'NumInterfaceNodes',din);
            end
            if ~isempty(dnc),
                dim_cons=netcdf.defDim(obj.Nid,'NumConstraints',dnc);
            end
            if deig, %% if deig is zero, don't write it out.
                dim_eig=netcdf.defDim(obj.Nid,'NumEig',deig);
            end
            
            if ~isempty(otm),
                dim_otm=netcdf.defDim(obj.Nid,'OTM_col',size(otm,1));
            end
            if ~isempty(outmap),
                dim_out=netcdf.defDim(obj.Nid,'NumNode_out',length(outmap));
            end
            if ~isempty(otme),
                dim_otme=netcdf.defDim(obj.Nid,'OTME_col',size(otme,1));
            end
            if ~isempty(outemap),
                dim_oute=netcdf.defDim(obj.Nid,'NumElem_out',length(outemap));
            end
            
            [ju,type]=obj.setType(1);
            if ~isempty(Kr),
                vark=netcdf.defVar(obj.Nid,'Kr',type,[dim_dof,dim_dof]);
            end
            if ~isempty(Cr),
                varcr=netcdf.defVar(obj.Nid,'Cr',type,[dim_dof,dim_dof]);
            end
            if ~isempty(Mr),
                varm=netcdf.defVar(obj.Nid,'Mr',type,[dim_dof,dim_dof]);
            end
            if ~isempty(cbmap),
                varc=netcdf.defVar(obj.Nid,'cbmap','int',[dim_2,dim_dof]);
            end
            if ~isempty(otm),
                varo=netcdf.defVar(obj.Nid,'OTM',type,[dim_dof,dim_otm]);
            end
            if ~isempty(outmap),
                varout=netcdf.defVar(obj.Nid,'OutMap','int',dim_out);
            end
            if ~isempty(otme),
                varoe=netcdf.defVar(obj.Nid,'OTME',type,[dim_dof,dim_otme]);
            end
            if ~isempty(outemap),
                varoute=netcdf.defVar(obj.Nid,'OutElemMap','int',dim_oute);
            end
            %%%%
            netcdf.endDef(obj.Nid);
            %%%%
            if ~isempty(Kr),
                netcdf.putVar(obj.Nid,vark,Kr');
            end
            if ~isempty(Cr),
                netcdf.putVar(obj.Nid,varcr,Cr');
            end
            if ~isempty(Mr),
                netcdf.putVar(obj.Nid,varm,Mr');
            end
            if ~isempty(cbmap),
                netcdf.putVar(obj.Nid,varc,cbmap');
            end
            if ~isempty(otm),
                netcdf.putVar(obj.Nid,varo,otm);
            end
            if ~isempty(outmap),
                netcdf.putVar(obj.Nid,varout,outmap);
            end
            if ~isempty(otme),
                netcdf.putVar(obj.Nid,varoe,otme);
            end
            if ~isempty(outemap),
                netcdf.putVar(obj.Nid,varoute,outemap);
            end
        end
        %%%%%
    end
    methods (Access=private)
        function [flag,varid]=isVar(obj,varname)
            flag=0;
            varid=-1;
            for i=1:obj.NumVariables,
                name=netcdf.inqVar(obj.Nid,i-1);
                if strcmpi(name,varname),
                    flag=1;
                    varid=i-1;
                    return
                end
            end
        end
        function [flag,attid]=isAtt(obj,varname)
            flag=0;
            attid=-1;
            for i=1:obj.NumGlobalAttributes,
                name=netcdf.inqAttName(obj.Nid,netcdf.getConstant('NC_GLOBAL'),i-1);
                if strcmpi(name,varname),
                    flag=1;
                    attid=i-1;
                    return
                end
            end
        end
        function [dimlen,dimid]=isDim(obj,varname)
            for i=1:obj.NumDimensions,
                [name,dimlen]=netcdf.inqDim(obj.Nid,i-1);
                if strcmpi(name,varname),
                    dimid=netcdf.inqDimID(obj.Nid,name);
                    return
                end
            end
           dimid=[];
           dimlen=[];
        end
        function [data,type]=setType(obj,data)
            switch obj.IO_WordSize,
                case 0
                    data=single(data);
                    type='float';
                case 4
                    data=single(data);
                    type='float';
                case 8
                    data=double(data);
                    type='double';
                otherwise
                    error('Unknown wordsize')
            end
        end
        function obj=updateProperties(obj)
            [obj.NumDimensions,obj.NumVariables,obj.NumGlobalAttributes,obj.UnlimitedDim]=netcdf.inq(obj.Nid);
        end
    end
    methods (Static,Access=private)
        function datanames=nem_datanames
            datanames={'el_blk_ids_global','el_blk_cnt_global','ns_ids_global','ns_node_cnt_global', ...
                'ns_df_cnt_global','ss_ids_global','ss_side_cnt_global','ss_df_cnt_global','nem_ftype', ...
                'int_n_stat','bor_n_stat','ext_n_stat','int_e_stat','bor_e_stat','elem_mapi','elem_mapb', ...
                'node_mapi','node_mapb','n_comm_ids','n_comm_stat','e_comm_ids','e_comm_stat','n_comm_data_idx', ...
                'n_comm_nids','n_comm_proc','e_comm_data_idx','e_comm_eids','e_comm_proc','e_comm_sids', ...
                [],[],[],[],[]};
        end
        function dimnames=nem_dimnames
            dimnames={'num_el_blk_global','num_el_blk_global','num_ns_global','num_ns_global', ...
                'num_ns_global','num_ss_global','num_ss_global','num_ss_global','','num_procs_file', ...
                'num_procs_file','num_procs_file','num_procs_file','num_procs_file','num_int_elem', ...
                'num_bor_elem','num_int_node','num_bor_node','num_n_cmaps','num_n_cmaps','num_e_cmaps', ...
                'num_e_cmaps','num_n_cmaps','ncnt_cmap','ncnt_cmap','num_e_cmaps','ecnt_cmap', ...
                'ecnt_cmap','ecnt_cmap','num_nodes_global','num_elems_global','num_ns_global', ...
                'num_ss_global','num_processors'};
        end
        function node=side_nodes(elemtype,elemcon,side,dim)
            % this is from exgssn.c in the exodus library
            % First convert element to blocknumber
            % then identify an element type with each element
            %
            switch lower(elemtype(1:3))
                case 'cir'
                    node=1;
                case 'sph'
                    node=1;
                case {'tru','bea'}
                    % truss and beams
                    node=[1 2 3];
                    if length(elemcon)==2,
                        node=node(1:2);
                    end
                case 'tri'
                    % triangle */
                    if dim==2,
                        switch side
                            case 1
                                node=[1,2,4];
                            case 2
                                node=[2,3,5];
                            case 3
                                node=[3,1,6];
                            otherwise
                                error('Error: Invalid %d face number %d',elemtype,side);
                        end
                        if length(elemcon)==3,
                            node=node(1:2);
                        end
                    else
                        % triangle 3d */ modified for 3 or 6 node only, could add middle node if necessary
                        switch side
                            case 1
                                node=[1,2,3,4,5,6]; %node=[1,2,3,4,5,6,7];
                            case 2
                                node=[3,2,1,6,5,4]; %node=[3,2,1,6,5,4,7];
                            case 3
                                node=[1,2,4];
                            case 4
                                node=[2,3,5];
                            case 5
                                node=[3,1,6];
                            otherwise
                                error('Error: Invalid %d face number %d',elemtype,side);
                        end
                        if length(elemcon)==3,
                            if side<3, % element face
                                node=node(1:3);
                            else
                                node=node(1:2);
                            end
                        end
                    end
                case 'qua'
                    % quad */
                    switch side
                        case 1
                            node=[1,2,5];
                        case 2
                            node=[2,3,6];
                        case 3
                            node=[3,4,7];
                        case 4
                            node=[4,1,8];
                        otherwise
                            error('Error: Invalid %d face number %d',elemtype,side);
                    end
                    if length(elemcon)==4,
                        node=node(1:2);
                    end
                case 'she'
                    % shell */
                    switch side
                        case 1
                            node=[1,2,3,4,5,6,7,8];
                        case 2
                            node=[1,4,3,2,8,7,6,5];
                        case 3
                            node=[1,2,5];
                        case 4
                            node=[2,3,6];
                        case 5
                            node=[3,4,7];
                        case 6
                            node=[4,1,8];
                        otherwise
                            error('Error: Invalid %d face number %d',elemtype,side);
                    end
                    if length(elemcon)==3,
                        if side<3, % element face
                            node=node(1:4);
                        else
                            node=node(1:2);
                        end
                    end
                case 'tet'
                    % tetra */
                    switch side
                        case 1
                            node=[1,2,4,5,9,8];
                        case 2
                            node=[2,3,4,6,10,9];
                        case 3
                            node=[1,4,3,8,10,7];
                        case 4
                            node=[1,3,2,7,6,5];
                        otherwise
                            error('Error: Invalid %d face number %d',elemtype,side);
                    end
                    if length(elemcon)==4,
                        node=node(1:3);
                    elseif length(elemcon)==8,
                        node=node(1:4);
                    end
                case 'wed'
                    % wedge */
                    switch side
                        case 1
                            node=[1,2,5,4,7,11,13,10];
                        case 2
                            node=[2,3,6,5,8,12,14,11];
                        case 3
                            node=[1,4,6,3,10,15,12,9];
                        case 4
                            node=[1,3,2,9,8,7];  % 3 node side
                        case 5
                            node=[4,5,6,13,14,15]; % 3 node side
                        otherwise
                            error('Error: Invalid %d face number %d',elemtype,side);
                    end
                    if length(elemcon)==6,
                        if side<4,
                            node=node(1:4);
                        else
                            node=node(1:3);
                        end
                    end
                case 'hex'
                    % hex */
                    switch side
                        case 1
                            node=[1,2,6,5,9,14,17,13,26];
                        case 2
                            node=[2,3,7,6,10,15,18,14,25];
                        case 3
                            node=[3,4,8,7,11,16,19,15,27];
                        case 4
                            node=[1,5,8,4,13,20,16,12,24];
                        case 5
                            node=[1,4,3,2,12,11,10,9,22];
                        case 6
                            node=[5,6,7,8,17,18,19,20,23];
                        otherwise
                            error('Error: Invalid %d face number %d',elemtype,side);
                    end
                    if length(elemcon)==8,
                        node=node(1:4);
                    elseif length(elemcon)==20,
                        node=node(1:8);
                    end
                case 'pyr'
                    % pyramid */
                    switch side
                        case 1
                            node=[1,2,5,6,11,10];
                        case 2
                            node=[2,3,5,7,12,11];
                        case 3
                            node=[3,4,5,8,13,12];
                        case 4
                            node=[1,5,4,10,13,9];
                        case 5
                            node=[1,4,3,2,9,8,7,6];
                        otherwise
                            error('Error: Invalid %d face number %d',elemtype,side);
                    end
                    if length(elemcon)==5
                        if side<5
                            node=node(1:3);
                        else
                            node=node(1:4);
                        end
                    end
                otherwise
                    error('Unknown element, %s',elemtype);
            end
        end
        function arr=mk_str_array(str_cell,len)
            n=length(str_cell);
            arr=repmat(char(0),[len,n]);
            for i=1:n,
                if ~isempty(str_cell{i})
                    if len < length(str_cell{i})
                        arr(:,i) = str_cell{i}(1:len)';
                    else
                        arr(1:length(str_cell{i}),i)=str_cell{i}';
                    end
                end
            end
        end  
    end
    
    
end