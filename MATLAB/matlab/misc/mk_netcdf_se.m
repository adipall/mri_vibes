function mk_netcdf_se(filename,Mass,Stiff,Damp,cbmap,grid_interface,...
            OTM, OtmMap, grid_otm)
% generate a netcdf formatted text file of a superelement.
% the file can be converted to binary form using ncgen
%
%mk_netcdf_se(filename,Mass,Stiff,Damp,cbmap,grid_interface,...
%            OTM, OtmMap, grid_otm)
%
% filename - the output file name (text)
% cbmap    - see below
% Mass     - required mass matrix
% Stiff    - required stiffness matrix
% Damp     - OPTIONAL damping matrix
% OTM      - OPTIONAL output transfer matrix
% precision- OPTIONAL number of digits in the output
% grid_interface - number of interface nodes X 4
%                  column 1 = grid id
%                  column 2 = x coordinate of interface grid
%                  column 3 = y coordinate of interface grid
%                  column 4 = z coordinate of interface grid
% grid_otm - ditto for output transfer matrix
%
% The CBMAP is a two column map from dofs to the grid and cid. The
% first column is the grid id. The second is the cid (1:6). The map
% is 0 0 for generalized dofs.
%
% The OtmMap is also two columns, and maps dofs in the OTM to grids/cid pairs
 
precision=15;

% standard exodus file definitions
len_string = 33;
len_line = 81;
num_dim = 3;

% output a single super element
num_elem = 1;
num_el_blk = 1;
num_el_in_blk1 = 1;
num_nod_per_el1 = 1; % may change below

% model dimensions
grids=[grid_interface' grid_otm']';
num_nodes=size(grids,1);
NumInterfaceNodes=size(grid_interface,1);
if ( NumInterfaceNodes>0 )
  num_nod_per_el1 = NumInterfaceNodes;
end
NumNode_out=size(grid_otm,1);
OTM_col=size(OTM,1);
% NumEig   - number of eigenvalues (fixed interface modes) in the matrices
NumConstraints=size(cbmap,1) - sum(cbmap(:,1)==0);

NumCols = length(Mass);
NumDof = length(Stiff);
NumEig = NumDof - NumConstraints;
doDamp=1;
if size(Damp,1) ~= NumDof
    doDamp=0;
end

%check a few dimensions
if size(cbmap,2)~=2
  error('cbmap must have two rows.');
end
if ( size(Mass,1) ~= size(Mass,2) || size(Stiff,1) ~= size(Stiff,2) )
  error('mass and stiffness matrices must be square.');
end
if NumCols ~=NumDof
    error('mass and stiffness dimensions are inconsistent %d != %d',NumCols,NumDof);
end
if OTM_col>0 && size(OtmMap,1)~= OTM_col
    error('OTM and OtmMap are inconsistent.\n');
end
if NumConstraints<=0
   error('no interface.\n');
end

% try to open the file and start the fun
fid=fopen(filename,'w');
if fid==-1
    error('unable to open file: "%s"',filename);
end

fprintf(fid,'netcdf test {\n');
disp(['writing text dumped netcdf data to ' filename]);

%Write Dimensions data:
fprintf(fid,'dimensions:\n');
fprintf(fid,['\tlen_string = ',num2str(len_string),' ;\n']);
fprintf(fid,['\tlen_line = ',num2str(len_line),' ;\n']);
fprintf(fid,'\ttime_step = UNLIMITED ;\n');
fprintf(fid,['\tnum_dim = ',num2str(num_dim),' ;\n']);
if ( num_nodes>0 )
  fprintf(fid,['\tnum_nodes = ',num2str(num_nodes),' ;\n']);
end
fprintf(fid,['\tnum_elem = ',num2str(num_elem),' ;\n']);
fprintf(fid,['\tnum_el_blk = ',num2str(num_el_blk),' ;\n']);
% fprintf(fid,['\tnum_node_sets = ',num2str(num_node_sets),' ;\n']);
fprintf(fid,['\tnum_el_in_blk1 = ',num2str(num_el_in_blk1),' ;\n']);
fprintf(fid,['\tnum_nod_per_el1 = ',num2str(num_nod_per_el1),' ;\n']);
fprintf(fid,['\tNumDof = ',num2str(NumDof),' ;\n']);
fprintf(fid,'\ttwo = 2 ;\n');
fprintf(fid,['\tNumConstraints = ',num2str(NumConstraints),' ;\n']);
fprintf(fid,['\tNumInterfaceNodes = ',num2str(NumInterfaceNodes),' ;\n']);
fprintf(fid,['\tNumEig = ',num2str(NumEig),' ;\n']);
fprintf(fid,['\tNumNode_out = ',num2str(NumNode_out),' ;\n']);
if ( OTM_col>0 )
  fprintf(fid,['\tOTM_col = ',num2str(OTM_col),' ;\n']);
end


%Write Variables data:
fprintf(fid,'variables:\n');
fprintf(fid,'\tdouble time_whole(time_step) ;\n');
fprintf(fid,'\tint eb_status(num_el_blk) ;\n');
fprintf(fid,'\tint eb_prop1(num_el_blk) ;\n');
fprintf(fid,'\t        eb_prop1:name = "ID" ;\n');
fprintf(fid,'\tdouble coord(num_dim,num_nodes) ;\n');
fprintf(fid,'\tchar coor_names(num_dim,len_string) ;\n');
fprintf(fid,'\tint connect1(num_el_in_blk1,num_nod_per_el1) ;\n');
fprintf(fid,'\t\t connect1:elem_type = "SUPER" ;\n');
fprintf(fid,'\tdouble Kr(NumDof,NumDof) ;\n');
fprintf(fid,'\tdouble Mr(NumDof,NumDof) ;\n');
if doDamp==1
   fprintf(fid,'\tdouble Cr(NumDof,NumDof) ;\n');
end
if ( OTM_col>0 )
  fprintf(fid,'\tdouble OTM(OTM_col,NumDof) ;\n');
end
fprintf(fid,'\tint cbmap(NumDof,two) ;\n');
% fprintf(fid,'\tint node_num_map(num_nodes) ;\n');

%Write Global attributes:
fprintf(fid,'\n// global attributes:\n');
fprintf(fid,'\t\t:api_version = 4.21f ;\n');
fprintf(fid,'\t\t:version = 3.01f ;\n');
fprintf(fid,'\t\t:floating_point_word_size = 8 ;\n');
fprintf(fid,'\t\t:file_size = 0 ;\n');
fprintf(fid,'\t\t:title = "cbr model" ;\n');


%Begin Data:
fprintf(fid,'data:\n');

fprintf(fid,' eb_status = 1 ;\n\n');
fprintf(fid,' eb_prop1 = 1 ;\n\n');
fprintf(fid,'\n coord = \n');
% coord = grids(:,2:4);
writearray(grids(:,2:4),len_line,precision,fid);

fprintf(fid,'\n coor_names =\n');
fprintf(fid,'  "coordx",\n');
fprintf(fid,'  "coordy",\n');
fprintf(fid,'  "coordz" ;\n');

fprintf(fid,'\n connect1 = \n');
writearray(1,len_line,precision,fid);

fprintf(fid,'\n Kr = \n');
writearray(Stiff,len_line,precision,fid);

fprintf(fid,'\n Mr = \n');
writearray(Mass,len_line,precision,fid);

if doDamp == 1
    fprintf(fid,'\n Cr = \n');
    writearray(Damp,len_line,precision,fid);
end
    
if  OTM_col>0 
   fprintf(fid,'\n OTM = \n');
   writearray(OTM,len_line,precision,fid);
end

fprintf(fid,'\n cbmap = \n');
writearray(cbmap',len_line,precision,fid);

fprintf(fid,'}\n');

fclose(fid);

disp(['Text dump to ' filename ' completed.']);
end % end of the mk_netcdf_se function

%--------------------------------------------------------------------------
% function only visible here
% writes data out, but wraps the data as needed

function writearray(array,len_line,precision,fid)
%
%

s=size(array);
len_array = s(1)*s(2);
row = 1;
while row<(len_array)
    [I,J]=ind2sub(s,row);
    Qtest=[' ',num2str(array(I,J),precision)];
    [I,J]=ind2sub(s,row+1);
    next = num2str(array(I,J),precision);
    while ((length(Qtest)+length(next))<(len_line-4))
        Qtest=[Qtest,', '];
        row=row+1;
        [I,J]=ind2sub(s,row);
        Q=[Qtest,num2str(array(I,J),precision)];
        Qtest = Q;
        if row+1<len_array
            [I,J]=ind2sub(s,row+1);
            next = num2str(array(I,J),precision);
        else
            next = ones(1,len_line);
        end
    end
    if row<len_array
        fprintf(fid,[Q,', \n']);
    elseif row==len_array
        fprintf(fid,[Q,' ;\n']);
    end
    row=row+1;
end
if row==len_array
    [I,J]=ind2sub(s,row);
    fprintf(fid,[' ',num2str(array(I,J),precision),' ;\n']);
end

end % end of the writearray function
