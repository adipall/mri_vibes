function [err]=copy_netcdf(filename,destfilename)
    %
    % This function will copy a netcdf file to a new netcdf file.
    %  This could be done with a system call but this will do it internal to matlab
    %
    %  err=copy_netcdf(filename, destination_filename);
    %
    %
    
    % Author Todd Simmermacher
    % Tel: 844-0096
    % Last Updated: 04/13/2010
    % Sandia National Laboratories
    %
    
    err=0;
    fid=netcdf.open(filename,'nowrite');
    [ndims,nvars,ngatts,unlimitDim]=netcdf.inq(fid);
    b64=netcdf.getConstant('64bit_offset');
    nc=netcdf.getConstant('noclobber');
    %share=netcdf.getConstant('nc_share');
    %%
    fid_new=netcdf.create(destfilename,bitor(b64,nc));
    %
    %  Deal with Global Attributes first
    gbl=netcdf.getConstant('Global');
   
    for i=0:ngatts-1,
        attname=netcdf.inqAttName(fid,gbl,i);
        netcdf.copyAtt(fid,gbl,attname,fid_new,gbl);    
    end
    %
    % Now Deal with Dimensions first
    for i=0:ndims-1
        [dimname,dimlen]=netcdf.inqDim(fid,i);
        netcdf.defDim(fid_new,dimname,dimlen);
    end
    %% Now define variables.  Two loops so we only have to go out of define mode once
    for i=0:nvars-1,
        [varname,xtype,dimids,natts(i+1)] = netcdf.inqVar(fid,i);
        varid = netcdf.defVar(fid_new,varname,xtype,dimids);
        if natts(i+1)>0,
            for j=0:natts(i+1)-1,
                attname=netcdf.inqAttName(fid,i,j);
                netcdf.copyAtt(fid,i,attname,fid_new,i)
            end
            
        end
    end
    netcdf.endDef(fid_new)
    for i=0:nvars-1,
        data = netcdf.getVar(fid,i);
        if ~isempty(data),
            netcdf.putVar(fid_new,i,data);
        end
    end
    
    
    netcdf.close(fid);
    netcdf.close(fid_new);
end
