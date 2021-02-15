function [fmap,result,err]=exo_map(code,map,ffrom,fto)
    result='';
    err='';
    switch lower(code)
        case 'mapvar'
            [fmap,result,err]=mapvar(map,ffrom,fto);
        case 'encore'
            [fmap,result,err]=encore(map,ffrom,fto);
            
    end
    
    
end

function [fmap,result,err]=mapvar(map,ffrom,fto)
    
    if ispc,
        fmap=-1;
        error('Mapvar is not supported in Windows')
    end
    rnd=randi(1e6,1);
    if ischar(ffrom)
        system(sprintf('cp %s EXO%d.e',ffrom,rnd));
    else
        exo_wr(sprintf('EXO%d.e',rnd),ffrom);
    end
    
    if ischar(ffrom)
        system(sprintf('cp %s EXO%d.g',fto,rnd));
    else
        exo_wr(sprintf('EXO%d.g',rnd),fto);
    end
    
    exfile=sprintf('MAPVAR%d.inp',rnd);
    
    
    fid=fopen(exfile,'w');
    
    for i=1:length(map),
        if isfield(map,'blocksfrom') && ~ischar(map(i).blocksfrom),
            bf=sprintf('%20d',map(i).blocksfrom);
        else
            bf='all';
        end
        if isfield(map,'blocksto') && ~ischar(map(i).blocksto),
            bt=sprintf('%20d',map(i).blocksto);
        else
            bt='all';
        end
        if isfield(map,'scheme'),
            sc=map(i).scheme;
        else
            sc=1;
        end
        
        if isfield(map,'deform')
            de=map(i).deform;
        else
            de=1;
        end
        fprintf(fid,'deform %d\n',de);
        
        if isfield(map,'timestep') && ~ischar(map(i).timestep),
            ts=sprintf('%20d',map(i).timestep);
        else
            ts='all';
        end
        fprintf(fid,'times %s\n',ts);
        
        fprintf(fid,'map %s to %s scheme %d\n',bf,bt,sc);
        fprintf(fid,'run\n');
        fprintf(fid,'exit\n');
        fclose(fid);
    end
    [result,err]=system(sprintf('mapvar EXO%d < MAPVAR%d.inp',rnd,rnd));
    
    if ~ischar(fto),
        fmap=exo_rd(sprintf('EXO%d.int',rnd));
    else
        fmap=sprintf('mapped_%s',fto);
        if exist(fmap,'file');
            delete(fmap);
        end
        system(sprintf('cp EXO%d.int %s',rnd,fmap));
    end
    delete(sprintf('EXO%d.*',rnd));
    delete(sprintf('MAPVAR%d.inp',rnd));
    
end


function [fmap,result,err]=encore(map,ffrom,fto)
    
    if ispc,
        fmap=-1;
        error('Encore is not supported in Windows')
    end
    
    fid=fopen(exfile,'w');
    rnd=randi(1e6,1);
    fprintf(fid,'Begin Sierra Encore\n\n');
    if ischar(ffrom),
        from=ffrom;
        filefrom=ffrom;
        obj=ExodusIO(filefrom,'rw');
        en=obj.rd_ElemVarsNames;
        blkidsfrom=obj.rd_BlockIDs;
        volflag=0;
        if strmatch(en,'VOL');
            volflag=1;
        end
        close(obj);
    else
        from=ffrom.Title;
        filefrom=sprintf('EXO%d.e',rnd);
        %% write the file
    end
    if ischar(fto),
        to=fto;
    else
        to=fto.Title;
        fileto=sprintf('EXO%d.e',rnd);
        %% write the file
    end
    fprintf(fid,'    TITLE Mapping from %s to %s\n\n',from,to);
    fprintf(fid,'    BEGIN Postprocessor Output Control pp_out\n');
    fprintf(fid,'        Output To Console\n');
    fprintf(fid,'        Write To File encoreinfo.txt\n');
    fprintf(fid,'    END\n\n');
    fprintf(fid,'    BEGIN Finite Element Model meshFrom\n');
    fprintf(fid,'        Database Name = %s\n',filefrom);
    fprintf(fid,'    END\n\n');
    fprintf(fid,'    BEGIN Finite Element Model meshTo\n');
    fprintf(fid,'        Database Name = %s\n',fileto);
    fprintf(fid,'    END\n\n');
    
    if volflag,
        %% need to be able to calculate VOL*EDEP/realvol
    end
    
    fprintf(fid,'    BEGIN Field Function eFrom\n');
    fprintf(fid,'        Use Element Field EDEP\n');
    fprintf(fid,'    END\n\n');
    
    fprintf(fid,'    BEGIN Encore Procedure trans\n');
    fprintf(fid,'        Use System main\n');
    fprintf(fid,'        BEGIN System main\n');
    fprintf(fid,'            BEGIN Transient etrans\n');
    fprintf(fid,'                Advance regionFrom\n');
    fprintf(fid,'                Transfer fromto\n');
    fprintf(fid,'                Advance regionTo\n');
    fprintf(fid,'            END\n');
    fprintf(fid,'        Simulation Start Time = 0\n');
    fprintf(fid,'        Simulation Termination Time = 1\n');
    fprintf(fid,'        Simulation Max Global Iterations = 100\n');
    fprintf(fid,'        END\n');
    fprintf(fid,'    END\n\n');
    
    
    for i=1:length(map),
        if isfield(map,'blocksfrom') && ~ischar(map(i).blocksfrom),
            bf=sprintf('%20d',map(i).blocksfrom);
        else
            bf='all';
        end
        if isfield(map,'blocksto') && ~ischar(map(i).blocksto),
            bt=sprintf('%20d',map(i).blocksto);
        else
            bt='all';
        end
        if isfield(map,'scheme'),
            sc=map(i).scheme;
        else
            sc=1;
        end
        
        if isfield(map,'deform')
            de=map(i).deform;
        else
            de=1;
        end
        fprintf(fid,'deform %d\n',de);
        
        if isfield(map,'timestep') && ~ischar(map(i).timestep),
            ts=sprintf('%20d',map(i).timestep);
        else
            ts='all';
        end
        fprintf(fid,'times %s\n',ts);
        
        fprintf(fid,'map %s to %s scheme %d\n',bf,bt,sc);
        fprintf(fid,'run\n');
        fprintf(fid,'exit\n');
        fclose(fid);
    end
    [result,err]=system(sprintf('mapvar EXO%d < MAPVAR%d.inp',rnd,rnd));
    
    if ~ischar(fto),
        fmap=exo_rd(sprintf('EXO%d.int',rnd));
    else
        fmap=sprintf('mapped_%s',fto);
        if exist(fmap,'file');
            delete(fmap);
        end
        system(sprintf('cp EXO%d.int %s',rnd,fmap));
    end
    delete(sprintf('EXO%d.*',rnd));
    delete(sprintf('MAPVAR%d.inp',rnd));
    
end